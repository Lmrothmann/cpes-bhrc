################################################################################
# 10_final_tables_corrected.R
# Psychosocial exposome (cPES) study — final, publication-ready summary tables.
#
# The model-fitting scripts (02_longitudinal, 03_cortical_thickness, ...) produce
# per-term estimates (beta, SE, t, p) that are correct. This script applies the
# CORRECT multiple-comparison families and the interpretable effect metric used
# in the manuscript:
#   (i)  cortical FDR = Benjamini-Hochberg WITHIN hemisphere over the 34
#        Desikan-Killiany parcels only (global volumes eTIV / BrainSeg / mean
#        thickness are NOT members of the family);
#   (ii) cognition FDR = BH over the FOUR exposome-score coefficients only
#        (intercepts / age / covariate rows excluded);
#   (iii) developmental-years equivalent, with 95% CI (delta method).
#
# Inputs: the per-term CSVs written by scripts 03 (cortex) and 02 (cognition).
# Reproduces: eTable S5 (cortex q34), eTable S6 (dev-years), eTable S7 (cognition).
# R 4.4.1. No random seed needed (deterministic recomputation from saved estimates).
################################################################################

suppressWarnings(suppressMessages({library(readr); library(dplyr)}))

## ---- configure paths (edit ROOT to your local output tree) -------------------
ROOT   <- Sys.getenv("cPES_OUT", "FINAL_WITH_SUBSTANCE")
CORTEX_MAIN <- file.path(ROOT, "03_cPES_BRAIN/LMR_EXTRA_THICKNESS_cPES_FINAL/01_MAIN_MODEL/02_MAIN_all_terms_FDR_LH_RH.csv")
CORTEX_ADJ  <- file.path(ROOT, "03b_cPES_BRAIN_ADJ_CBCL/LMR_EXTRA_THICKNESS_cPES_FINAL/01_MAIN_MODEL/02_MAIN_all_terms_FDR_LH_RH.csv")
COG_TIDY    <- file.path(ROOT, "02_LONGITUDINAL_COGNITION/all_models_tidy_corrected.csv")
OUT         <- file.path(ROOT, "_CORRECTED_2026-07"); dir.create(OUT, showWarnings = FALSE)

NON_CORTICAL <- c("eTIV", "BrainSegVolNotVent")          # excluded from the parcel family
is_mean      <- function(x) grepl("MeanThickness", x)

## ---- (i) CORTEX: age x cPES, FDR over 34 parcels within hemisphere -----------
cortex_fdr <- function(path, label){
  d <- read_csv(path, show_col_types = FALSE)
  d <- d[d$term == "age_wave:cPES_z" & !(d$roi %in% NON_CORTICAL) & !is_mean(d$roi), ]
  d <- do.call(rbind, lapply(split(d, d$table), function(h){
        h$pFDR_34 <- p.adjust(h$p.value, "BH"); h[order(h$p.value), ] }))
  d$analysis <- label
  d$tier <- ifelse(d$pFDR_34 < .05, "confirmatory",
             ifelse(d$pFDR_34 < .10, "exploratory", "ns"))
  d[, c("analysis","table","roi","estimate","std.error","statistic",
        "p.value","pFDR_34","tier","n_sub","n_obs")]
}
cortex <- rbind(cortex_fdr(CORTEX_MAIN, "MAIN_unadjusted"),
                cortex_fdr(CORTEX_ADJ,  "MAIN_adjBaselineCBCL"))
write_csv(cortex, file.path(OUT, "02_CORTEX_FDR_over34parcels.csv"))

## ---- (ii) DEVELOPMENTAL-YEARS EQUIVALENT with 95% CI (delta method) ----------
# additional thinning per interquartile increase in cPES, over a 5-year window,
# divided by the region's normative annual thinning rate.
DELTA_cPES <- 1.4724632884458098   # IQR of cPES (z): Q3 - Q1
WINDOW     <- 5                     # years
K          <- DELTA_cPES * WINDOW   # = 7.362

devyears <- function(path, label){
  d <- read_csv(path, show_col_types = FALSE)
  d <- d[!(d$roi %in% NON_CORTICAL) & !is_mean(d$roi), ]
  ax <- d[d$term == "age_wave:cPES_z", ]; ag <- d[d$term == "age_wave", ]
  m <- merge(ax[, c("table","roi","estimate","std.error","p.value")],
             ag[, c("table","roi","estimate")], by = c("table","roi"),
             suffixes = c("_ageXcPES","_age"))
  m <- do.call(rbind, lapply(split(m, m$table), function(h){
        h$pFDR_34 <- p.adjust(h$p.value, "BH"); h }))
  m$dev_years <- (m$estimate_ageXcPES * K) / abs(m$estimate_age)
  # 95% CI: propagate SE of the interaction (age main-effect ~4% rel.SE, negligible)
  se_ratio <- abs(m$estimate_ageXcPES / m$estimate_age) * (m$std.error / abs(m$estimate_ageXcPES))
  m$dev_years_lo <- m$dev_years - 1.96 * K/abs(m$estimate_age) * m$std.error
  m$dev_years_hi <- m$dev_years + 1.96 * K/abs(m$estimate_age) * m$std.error
  m$analysis <- label
  m[order(m$table, m$p.value), ]
}
dy <- rbind(devyears(CORTEX_MAIN, "MAIN_unadjusted"),
            devyears(CORTEX_ADJ,  "MAIN_adjBaselineCBCL"))
write_csv(dy, file.path(OUT, "03_developmental_years_RECOMPUTED.csv"))

## ---- (iii) COGNITION: FDR over the four exposome coefficients -----------------
cg <- read_csv(COG_TIDY, show_col_types = FALSE)
cog_family <- function(term_keep, label){
  s <- cg[cg$term == term_keep, ]
  s$FDR_4tests <- p.adjust(s$p.value, "BH")
  s$model <- label; s$term_kept <- term_keep
  s[order(s$p.value), c("model","term_kept","Outcome","estimate","std.error",
                        "p.value","FDR_4tests")]
}
cognition <- rbind(cog_family("env_c", "UNADJ_main"),
                   cog_family("time:env_c", "UNADJ_timeXcpes"))
cognition$tier <- ifelse(cognition$FDR_4tests < .05, "confirmatory",
                   ifelse(cognition$FDR_4tests < .10, "exploratory", "ns"))
write_csv(cognition, file.path(OUT, "01_COGNITION_FDR_over4tests.csv"))

## ---- console summary ---------------------------------------------------------
cat("=== CORTEX age x cPES (right hemisphere, q over 34 parcels) ===\n")
print(cortex[cortex$table == "rh_aparc_thickness" & cortex$analysis == "MAIN_unadjusted",
             c("roi","estimate","p.value","pFDR_34","tier")][1:3, ], row.names = FALSE)
cat("\n=== DEVELOPMENTAL YEARS (unadjusted, 3 right-hemisphere confirmatory ROIs) ===\n")
key_rh <- c("rh_inferiorparietal_thickness","rh_supramarginal_thickness","rh_lateralorbitofrontal_thickness")
print(dy[dy$analysis == "MAIN_unadjusted" & dy$roi %in% key_rh,
         c("roi","dev_years","dev_years_lo","dev_years_hi")], row.names = FALSE)
cat("\n=== COGNITION env_c (FDR over 4 tests) ===\n")
print(cognition[cognition$model == "UNADJ_main", c("Outcome","estimate","FDR_4tests","tier")], row.names = FALSE)
cat("\nWrote 3 tables to:", OUT, "\n")
