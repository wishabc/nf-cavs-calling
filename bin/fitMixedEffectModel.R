library(lme4)
library(lmerTest)
library(variancePartition)
library(data.table)


weighted.var <- function(x, w, na.rm = FALSE) {
    if (missing(w)) {
    w <- rep.int(1, length(x))
    } else if (length(w) != length(x)) {
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

process_group <- function(current_data, vpcontrol) {
    # Calculate weights within the group
    current_data$w <- current_data$inverse_mse / mean(current_data$inverse_mse)

    full_variance <- weighted.var(current_data$es, current_data$w)

    placeholder_df <- expand.grid(
        group_id=unique(current_data$group_id),
        group_es=NA_real_,
        group_es_std=NA_real_,
        indiv_id_rand_var=NA_real_,
        chisq=NA_real_,
        chi_df=NA_real_,
        p_differential=NA_real_,
        var_indiv_id = NA_real_,
        var_group_id = NA_real_,
        var_residuals = NA_real_,
        es_var = NA_real_,
        samples_count = NA_integer_,
        full_variance = full_variance
    )

    tryCatch({
        # Fit the linear mixed-effects models
        full_model <- lmer(es ~ 0 + group_id + (1 | indiv_id), 
            data = current_data, weights=w, REML=FALSE,
            control=vpcontrol
        )

        reduced_model <- lmer(es ~ 1 + (1 | indiv_id),
            data = current_data, weights=w, REML=FALSE,
            control=vpcontrol
        )

        coef_df <- as.data.frame(summary(full_model)$coefficients)

        # Add group_id column based on the names of the coefficients
        coef_df$group_id <- gsub("group_id", "", rownames(coef_df))
        merged_data <- as.data.frame(merge(current_data, coef_df, by="group_id"))
        print(str(merged_data))
        merged_data$es_diff <- merged_data$es - merged_data$"Estimate"
        es_var <- aggregate(
            merged_data[c('es_diff', 'w')], 
            list(merged_data$group_id),
            FUN=function(x) {
                w_var = weighted.var(x$es_diff, x$w)
                count = length(x$es_diff)
                return(c(w_var, count))
            }
        )
        colnames(es_var) <- c("group_id", "es_var", "samples_count")
        coef_df <- merge(coef_df, es_var, by="group_id")

        # Extract variance and standard deviation of random effect for indiv_id
        random_effect_variance_indiv_id <- VarCorr(full_model)$indiv_id[1,1]

        # Extract residual variance
        residual_variance <- attr(lme4::VarCorr(full_model), "sc")^2

        # Add random effect variance, standard deviation, and residual variance to the dataframe
        coef_df$indiv_id_rand_var <- rep(random_effect_variance_indiv_id, nrow(coef_df))

        # ANOVA between full and reduced models
        anova_result <- anova(reduced_model, full_model)
        chisq <- anova_result$"Chisq"[2]
        df <- anova_result$"Df"[2]
        p_value <- anova_result$"Pr(>Chisq)"[2]

        # Add these values as new columns to coef_df, repeated for each row
        coef_df$chisq <- rep(chisq, nrow(coef_df))
        coef_df$chi_df <- rep(df, nrow(coef_df))
        coef_df$p_differential <- rep(p_value, nrow(coef_df))

        # Calculating variance
        fit <- full_model
        varComp <- lapply(lme4::VarCorr(fit), function(fit) attr(fit, "stddev")^2)
        varComp$Residuals <- sigma(fit)^2

        idx <- which(colnames(fit@pp$X) != "(Intercept)")
            # print(paste("Length of idx:", length(idx)))
            # if (length(idx) == 0) {
            #     print("No fixed effects other than intercept found.")
            # }
        fxeff <- sapply(idx, function(i) {
            fit@pp$X[, i] * fixef(fit)[i]
        })
        colnames(fxeff) <- colnames(fit@pp$X)[idx]

        N <- nrow(fxeff)
        varFixedTotal <- weighted.var(rowSums(fxeff), w) * (N - 1) / N
        varComp[['group_id']] <- varFixedTotal

        vc <- unlist(varComp)
        res <- vc / sum(vc)

        # remove ".(Intercept)" string
        names(res) <- gsub("\\.\\(Intercept\\)", "", names(res))

        coef_df$var_indiv_id <- res[['indiv_id']]
        coef_df$var_group_id <- res[['group_id']]
        coef_df$var_residuals <- res[['Residuals']]
        coef_df$full_variance <- full_variance

        return(coef_df)
    }, warning = function(w) {

        print(paste("Warning in model fit:", w$message))
        return(placeholder_df)
    }, error = function(e) {

        print(paste("Error in model fit:", e$message))
        return(placeholder_df)
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

# Split the data by the grouping columns
results <- snps_df[, process_group(.SD, vpcontrol), by='variant_id']  

# Write the final dataframe to a file
fwrite(results, file=outpath, sep = "\t", row.names = FALSE, quote = FALSE)