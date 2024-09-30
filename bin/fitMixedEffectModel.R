library(lme4)
library(lmerTest)
library(data.table)


weighted.var <- function(x, w, na.rm = FALSE) {
    if (missing(w)) {
        w <- rep.int(1, length(x))
    }
    
    if (length(w) != length(x)) {
        stop("x and w must have the same length")
    }

    if (min(w) < 0) {
        stop("there are negative weights")
    }

    if (is.integer(w)) {
        w <- as.numeric(w)
    }
    if (na.rm) {
        w <- w[obs.ind <- !is.na(x)]
        x <- x[obs.ind]
    }
    w <- w * length(w) / sum(w)
    return(sum(w * (x - weighted.mean(x, w))^2) / (sum(w) - 1))
}

vpcontrol <- lme4::lmerControl(
    calc.derivs = FALSE,
    check.rankX = "stop.deficient",
    check.conv.singular = lme4::.makeCC("ignore", tol = 1e-4)
)

process_group <- function(testedDF, vpcontrol) {
    # Calculate weights within the group
    testedDF$w <- testedDF$inverse_mse / mean(testedDF$inverse_mse)

    full_variance <- weighted.var(testedDF$es, testedDF$w)

    placeholderDF <- expand.grid(
        group_id=unique(testedDF$group_id),
        group_es=NA_real_,
        group_es_std=NA_real_,
        df=NA_real_,    
        "t value"=NA_real_,
        "Pr(>|t|)"=NA_real_,
        indivs_count = NA_integer_,
        samples_count = NA_integer_,
        indiv_id_rand_var=NA_real_,
        chisq=NA_real_,
        chi_df=NA_real_,
        p_differential=NA_real_,
        var_indiv_id = NA_real_,
        var_group_id = NA_real_,
        var_residuals = NA_real_,
        full_variance = full_variance
    )
    placeholderDF$group_id <- as.character(placeholderDF$group_id)

    tryCatch({
        # Fit the linear mixed-effects models
        fullModel <- lmer(es ~ 0 + group_id + (1 | indiv_id), 
            data = testedDF, weights=w, REML=FALSE,
            control=vpcontrol
        )

        reducedModel <- lmer(es ~ 1 + (1 | indiv_id),
            data = testedDF, weights=w, REML=FALSE,
            control=vpcontrol
        )

        coefficientDF <- setDT(as.data.frame(summary(fullModel)$coefficients))

        # Add group_id column based on the names of the coefficients
        coefficientDF$group_id <- gsub("group_id", "", rownames(summary(fullModel)$coefficients))
        merged_data <- merge(testedDF, coefficientDF, by="group_id")
        # merged_data$es_diff <- merged_data$es - merged_data$"Estimate"
        es_var <- merged_data[, .(indivs_count = uniqueN(indiv_id),
                                samples_count = .N),
                    by = .(group_id)]
        coefficientDF <- merge(coefficientDF, es_var, by="group_id")

        names(coefficientDF)[names(coefficientDF) == "Estimate"] <- "group_es"
        names(coefficientDF)[names(coefficientDF) == "Std. Error"] <- "group_es_std"

        # Extract variance and standard deviation of random effect for indiv_id
        fitStats <- lme4::VarCorr(fullModel)
        random_effect_variance_indiv_id <- fitStats$indiv_id[1,1]

        # Extract residual variance
        residual_variance <- attr(fitStats, "sc")^2

        # Add random effect variance, standard deviation, and residual variance to the dataframe
        coefficientDF$indiv_id_rand_var <- rep(random_effect_variance_indiv_id, nrow(coefficientDF))

        # LRT between full and reduced models
        lrtResult <- anova(reducedModel, fullModel)
        chisq <- lrtResult$"Chisq"[2]
        df <- lrtResult$"Df"[2]
        Pvalue <- lrtResult$"Pr(>Chisq)"[2]

        # Add these values as new columns to coefficientDF, repeated for each row
        coefficientDF$chisq <- rep(chisq, nrow(coefficientDF))
        coefficientDF$chi_df <- rep(df, nrow(coefficientDF))
        coefficientDF$p_differential <- rep(Pvalue, nrow(coefficientDF))

        # Calculating variance
        varComp <- lapply(fitStats, function(fullModel) attr(fullModel, "stddev")^2)
        varComp$Residuals <- sigma(fullModel)^2

        idx <- which(colnames(fullModel@pp$X) != "(Intercept)")
            # print(paste("Length of idx:", length(idx)))
            # if (length(idx) == 0) {
            #     print("No fixed effects other than intercept found.")
            # }
        fxeff <- sapply(idx, function(i) {
            fullModel@pp$X[, i] * fixef(fullModel)[i]
        })
        colnames(fxeff) <- colnames(fullModel@pp$X)[idx]

        N <- nrow(fxeff)
        varFixedTotal <- weighted.var(rowSums(fxeff), testedDF$w) * (N - 1) / N
        varComp[['group_id']] <- varFixedTotal

        vc <- unlist(varComp)
        res <- vc / sum(vc)

        # remove ".(Intercept)" string
        names(res) <- gsub("\\.\\(Intercept\\)", "", names(res))

        coefficientDF$var_indiv_id <- res[['indiv_id']]
        coefficientDF$var_group_id <- res[['group_id']]
        coefficientDF$var_residuals <- res[['Residuals']]
        coefficientDF$full_variance <- full_variance

        return(coefficientDF)
    }, warning = function(w) {

        print(paste("Warning in model fullModel:", w$message))
        return(placeholderDF)
    }, error = function(e) {

        print(paste("Error in model fullModel:", e$message))
        return(placeholderDF)
    })
}

args=(commandArgs(TRUE))
if (length(args)==0) {
  stop("No arguments supplied.")
} else {
  inpath <- args[1]
  outpath <- args[2]
}

snps_df <- fread(inpath)
snps_df$indiv_id <- as.factor(snps_df$indiv_id)
snps_df$group_id <- as.factor(snps_df$group_id)

grouping_columns <- c(colnames(snps_df)[1:6], 'variant_id')
# Split the data by the grouping columns
results <- snps_df[, process_group(.SD, vpcontrol), by=grouping_columns]  

# Write the final dataframe to a file
fwrite(results, file=outpath, sep = "\t", row.names = FALSE, quote = FALSE)