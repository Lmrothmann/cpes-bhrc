################################################################################
# 12_incremental_dR2_bootstrap_CI.R
#
# Produces the cross-validated incremental R2 (dR2) of adding the baseline cPES
# to a baseline-CBCL model at each follow-up wave, AND its 95% bootstrap CI.
#
# This is the follow-up-CI calculation: for each of 2,000 participant-level
# bootstrap resamples, the full 10-fold cross-validation is re-run for BOTH the
# baseline-CBCL model and the baseline-CBCL + cPES model, and the difference in
# out-of-fold R2 is recorded. The CI is the 2.5/97.5 percentiles of that
# distribution. (The baseline cPES enters as a fixed, previously learned
# out-of-fold predictor; the cross-validation re-estimates only the regression
# coefficients, not the cPES itself.)
#
# Placeholders: set <BHRC_INPUT_DIR> and ${cPES_OUT}. R 4.4.1, seed 241094.
################################################################################
suppressWarnings(suppressMessages({ library(haven); library(dplyr); library(readr); library(janitor); library(caret) }))
set.seed(241094)

PATH_BHRC <- "<BHRC_INPUT_DIR>/Grassi-Oliveira_167BHRC_2023_12_09.sav"
PATH_cPES <- "${cPES_OUT}/01_EXWAS_EN_cPES/POOLED/07_cpes_oof_predictions.csv"

B <- haven::read_sav(PATH_BHRC, encoding = "latin1") |> janitor::clean_names()
cPES <- readr::read_csv(PATH_cPES, show_col_types = FALSE) |>
  dplyr::transmute(ident = as.character(ident), cPES = as.numeric(cPES_score_OOF_z))

wv <- c("wave0_arm_1" = 0, "wave1_arm_1" = 1, "wave2_arm_1" = 2)
D <- B |> dplyr::filter(redcap_event_name %in% names(wv)) |>
  dplyr::mutate(wave = wv[redcap_event_name], ident = as.character(ident),
                cl_tot = as.numeric(cl_tot), age = as.numeric(age))

## outcome transform fitted once on all CBCL observations, as in the paper
allc <- D |> dplyr::filter(!is.na(cl_tot))
pre  <- caret::preProcess(allc["cl_tot"], method = c("YeoJohnson", "center", "scale"))
D$cl_tr <- NA_real_; ok <- !is.na(D$cl_tot)
D$cl_tr[ok] <- predict(pre, D[ok, "cl_tot", drop = FALSE])[[1]]

base <- D |> dplyr::filter(wave == 0) |>
  dplyr::transmute(ident, cbcl0_tr = cl_tr, age0 = age,
                   gender = factor(gender), site = factor(site_original)) |>
  dplyr::inner_join(cPES, by = "ident")

cv_r2 <- function(form, data, K = 10) {
  data <- data[stats::complete.cases(data[, all.vars(form)]), ]
  folds <- caret::createFolds(seq_len(nrow(data)), k = K)
  yhat <- rep(NA_real_, nrow(data))
  for (f in folds) { fit <- lm(form, data = data[-f, ]); yhat[f] <- predict(fit, newdata = data[f, ]) }
  y <- data[[all.vars(form)[1]]]; 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
}
f2 <- y ~ cbcl0_tr + age_fu + gender + site                 # baseline CBCL model
f3 <- y ~ cbcl0_tr + cPES + age_fu + gender + site          # + baseline cPES
R  <- 2000

for (w in c(1, 2)) {
  fu <- D |> dplyr::filter(wave == w, !is.na(cl_tr)) |> dplyr::transmute(ident, y = cl_tr, age_fu = age)
  dd <- base |> dplyr::inner_join(fu, by = "ident") |>
    dplyr::filter(!is.na(cbcl0_tr), !is.na(cPES), !is.na(age_fu), !is.na(gender), !is.na(site))
  point <- cv_r2(f3, dd) - cv_r2(f2, dd)
  bs <- numeric(R)
  for (b in 1:R) { idx <- sample(nrow(dd), replace = TRUE); db <- dd[idx, ]; bs[b] <- cv_r2(f3, db) - cv_r2(f2, db) }
  ci <- quantile(bs, c(.025, .975), names = FALSE)
  cat(sprintf("WAVE %d (n = %d): cross-validated dR2 = %+.3f  95%% CI [%+.3f, %+.3f]\n",
              w, nrow(dd), point, ci[1], ci[2]))
}
################################################################################
