################################################################################
# apply_cPES.R — apply the fitted cPES to new participants.
#
# Requires the complete application parameters produced by 01_build_cPES.R:
#   cPES_apply_parameters.csv  (one row per exposure, with columns:
#       variable, impute_median, yeojohnson_lambda, center, scale, weight)
#   cPES_intercept.txt         (the calibration intercept, in standardized-outcome units)
#
# For each exposure, the full pipeline is reproduced on new data:
#   1) median imputation      -> impute_median
#   2) Yeo-Johnson transform  -> yeojohnson_lambda
#   3) center and scale       -> (x - center) / scale
#   4) weighted sum + intercept
#
# The 52 non-zero weights plus the intercept and these preprocessing parameters
# constitute the complete model required to compute the cPES for a new sample.
################################################################################
suppressWarnings(suppressMessages({ library(readr); library(dplyr) }))

# Yeo-Johnson transform (vector x, scalar lambda); NA lambda = identity.
yeojohnson <- function(x, lambda) {
  if (is.na(lambda)) return(x)
  out <- numeric(length(x)); pos <- !is.na(x) & x >= 0; neg <- !is.na(x) & x < 0
  if (lambda != 0) out[pos] <- ((x[pos] + 1)^lambda - 1) / lambda else out[pos] <- log(x[pos] + 1)
  if (lambda != 2) out[neg] <- -(((-x[neg] + 1)^(2 - lambda) - 1) / (2 - lambda)) else out[neg] <- -log(-x[neg] + 1)
  out
}

# new_data: data.frame with one row per participant and the exposure columns
#           (column names matching 'variable' in the parameter file).
apply_cPES <- function(new_data, params_csv = "cPES_apply_parameters.csv", intercept) {
  P <- readr::read_csv(params_csv, show_col_types = FALSE)
  missing_cols <- setdiff(P$variable, names(new_data))
  if (length(missing_cols)) warning("missing exposure columns: ", paste(missing_cols, collapse = ", "))
  score <- rep(intercept, nrow(new_data))
  for (i in seq_len(nrow(P))) {
    if (!P$variable[i] %in% names(new_data)) next
    x <- as.numeric(new_data[[P$variable[i]]])
    x[is.na(x)] <- P$impute_median[i]                 # 1) median imputation
    x <- yeojohnson(x, P$yeojohnson_lambda[i])        # 2) Yeo-Johnson
    x <- (x - P$center[i]) / P$scale[i]               # 3) center + scale
    score <- score + P$weight[i] * x                  # 4) weighted sum
  }
  score
}

# --- Example -----------------------------------------------------------------
# intercept <- as.numeric(sub(".*: ", "", readLines("cPES_intercept.txt")[1]))
# cpes <- apply_cPES(new_data, "cPES_apply_parameters.csv", intercept)
################################################################################
