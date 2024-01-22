library(lme4)
library(lmerTest)
library(variancePartition)
library(dplyr)
library(tidyr)
library(purrr)

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

args=(commandArgs(TRUE))
if (length(args)==0) {
  stop("No arguments supplied.")
} else {
  inpath <- args[1]
  outpath <- args[2]
}

snps_df <- read.delim(inpath)
snps_df

# Define the starting columns for grouping
starting_columns_names <- names(snps_df)[1:6]

# Split the data by the grouping columns
split_data <- split(snps_df, snps_df[, starting_columns_names])

split_data <- snps_df %>%
  group_by(across(all_of(starting_columns_names))) %>%
  group_split()

# Process each group
results <- map(split_data, ~ {
  current_data <- .x

  # Calculate weights within the group
  w <- current_data$inverse_mse / mean(current_data$inverse_mse)

  # Fit the linear mixed-effects model
  full_model <- lme4::lmer(es ~ 0 + group_id + (1 | indiv_id), data = current_data, weights = w, REML=FALSE, control=vpcontrol)

  # Extract fixed effects coefficients and standard errors
  coefficients <- fixef(full_model)
  se <- summary(full_model)$coefficients[, "Std. Error"]
  names(coefficients) <- gsub("group_id", "", names(coefficients))

  # Create a dataframe for coefficients and standard errors
  coef_df <- data.frame(group_id = names(coefficients), group_es = coefficients, group_es_std = se)

  # Extract variance and standard deviation of random effect for indiv_id
  random_effect_variance_indiv_id <- VarCorr(full_model)$indiv_id[1,1]

  # Extract residual variance
  residual_variance <- attr(VarCorr(full_model), "sc")^2

  # Add random effect variance, standard deviation, and residual variance to the dataframe
  coef_df$indiv_id_rand_var <- rep(random_effect_variance_indiv_id, nrow(coef_df))

  # Fit the linear mixed-effects model
  reduced_model <- lme4::lmer(es ~ (1 | indiv_id), data = current_data, weights=w, REML=FALSE, control=vpcontrol)
  anova_result <- anova(reduced_model, full_model)

  chisq <- anova_result$"Chisq"[2]  
  df <- anova_result$"Df"[2]
  p_value <- anova_result$"Pr(>Chisq)"[2]

  # Add these values as new columns to coef_df, repeated for each row
  coef_df$chisq <- rep(chisq, nrow(coef_df))
  coef_df$chi_df <- rep(df, nrow(coef_df))
  coef_df$p_differential <- rep(p_value, nrow(coef_df))

  fit <- reduced_model
  #fit <- refit(fit, control = vpcontrol)
  # summary(fit)

  varComp <- lapply(lme4::VarCorr(fit), function(fit) attr(fit, "stddev")^2)
  varComp$Residuals <- sigma(fit)^2
  # varComp

  idx <- which(colnames(fit@pp$X) != "(Intercept)")
  fxeff <- sapply(idx, function(i) {
    fit@pp$X[, i] * lme4::fixef(fit)[i]
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

  # Combine the first row of original data with the coefficients dataframe
  combined_df <- cbind(current_data[1, starting_columns_names, drop = FALSE], coef_df)

  return(combined_df)
})

# Combine all the results into one dataframe
final_df <- bind_rows(results)
names(final_df)[1] <- '#chr'

# Write the final dataframe to a file
write.table(final_df, file = outpath, sep = "\t", row.names = FALSE, quote = FALSE)