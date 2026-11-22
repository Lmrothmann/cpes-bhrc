################################################################################
# 11_sample_Ns_and_roi_table.R
# Cumulative psychosocial exposome score (cPES) study - BHRC.
# Analytic sample sizes per wave (CBCL, each cognitive test, MRI) and the per-ROI table.
#
# Faithful, path-genericised copy of the analysis code used for the manuscript.
# CONFIGURE PATHS: set the placeholders below to your local tree before running.
#   <PROJECT_ROOT>     project directory
#   ${cPES_OUT}        output tree (e.g. <PROJECT_ROOT>/FINAL_WITH_SUBSTANCE)
#   <BHRC_INPUT_DIR>   controlled BHRC data (.sav) - not shareable
#   <COGNITION_DIR>    cognitive-score file
# R 4.4.1; seed 241094 where stochastic.
# NB: any FDR columns written here use the per-model pool; the MANUSCRIPT FDR
#     families and the developmental-years metric are applied in
#     10_final_tables_corrected.R (run last). See README.md.
################################################################################

################################################################################
# Reviewer #10/#11: analytic N per wave (CBCL, each cognitive outcome, MRI)
# Reviewer #23/#24: per-ROI table (beta/SE/CI/p/FDR) unadjusted vs adjusted
################################################################################
suppressWarnings(suppressMessages({library(haven); library(dplyr); library(readr); library(janitor); library(tidyr)}))

PATH_BHRC <- "<BHRC_INPUT_DIR>/Grassi-Oliveira_167BHRC_2023_12_09.sav"
PATH_cPES  <- "${cPES_OUT}/01_EXWAS_EN_cPES/POOLED/07_cpes_oof_predictions.csv"
PATH_COG  <- "<COGNITION_DIR>/Exec"
WSB       <- "${cPES_OUT}"
OUT       <- file.path(WSB, "_REVIEW_ANALYSES"); if(!dir.exists(OUT)) dir.create(OUT, recursive=TRUE)

cPES <- read_csv(PATH_cPES, show_col_types=FALSE) %>% transmute(ident=as.character(ident), cPES_score_OOF_z)
B <- read_sav(PATH_BHRC, encoding="latin1") %>% clean_names()
wv <- function(x) c("wave0_arm_1"=0,"wave1_arm_1"=1,"wave2_arm_1"=2)[x]

## ---- CBCL per wave ----
cb <- B %>% filter(redcap_event_name %in% names(c("wave0_arm_1"=0,"wave1_arm_1"=1,"wave2_arm_1"=2))) %>%
  transmute(ident=as.character(ident), wave=wv(redcap_event_name), cl_tot=as.numeric(cl_tot),
            age=as.numeric(age), gender, site=site_original) %>%
  inner_join(cPES, by="ident") %>% filter(!is.na(cl_tot), !is.na(age), !is.na(gender), !is.na(site))
cbcl_N <- cb %>% count(wave, name="N_CBCL")

## ---- Cognition per wave (per outcome) ----
Ex <- readRDS(PATH_COG) %>% clean_names() %>% transmute(ident=as.character(ident), wave=wv(redcap_event_name), wm_fl, ic_fl, tp_fl, ef_fl)
cog <- B %>% filter(redcap_event_name %in% names(c("wave0_arm_1"=0,"wave1_arm_1"=1,"wave2_arm_1"=2))) %>%
  transmute(ident=as.character(ident), wave=wv(redcap_event_name), age=as.numeric(age), gender, site=site_original) %>%
  inner_join(cPES, by="ident") %>% inner_join(Ex, by=c("ident","wave")) %>%
  filter(!is.na(age), !is.na(gender), !is.na(site))
cog_N <- cog %>% pivot_longer(c(wm_fl,ic_fl,tp_fl,ef_fl), names_to="outcome", values_to="v") %>%
  filter(!is.na(v)) %>% count(outcome, wave) %>% pivot_wider(names_from=wave, values_from=n, names_prefix="wave")

## ---- MRI per wave (from thickness QC) ----
qc <- read_csv(file.path(WSB,"03_cPES_BRAIN/LMR_EXTRA_THICKNESS_cPES_FINAL/00_QC/00_QC_ROI_SAMPLE_OVERVIEW.csv"), show_col_types=FALSE)
mri_N <- qc %>% filter(grepl("thickness", table)) %>%
  summarise(n_sub=max(n_sub), n_obs=max(n_obs), wave0=max(n_wave0), wave1=max(n_wave1), wave2=max(n_wave2))

cat("=== CBCL N per wave ===\n"); print(cbcl_N, row.names=FALSE)
cat("\n=== Cognition N per wave (per outcome) ===\n"); print(as.data.frame(cog_N), row.names=FALSE)
cat("\n=== MRI (thickness) N ===\n"); print(as.data.frame(mri_N), row.names=FALSE)
write_csv(cbcl_N, file.path(OUT,"10_11_CBCL_N_per_wave.csv"))
write_csv(cog_N, file.path(OUT,"10_11_cognition_N_per_wave.csv"))
write_csv(mri_N, file.path(OUT,"10_11_MRI_N.csv"))

## ---- Per-ROI table: unadjusted vs adjusted age×cPES ----
tg <- c("rh_lateralorbitofrontal_thickness","rh_inferiorparietal_thickness","rh_supramarginal_thickness")
grab <- function(path, tag){
  m <- read_csv(path, show_col_types=FALSE)
  m %>% filter(term=="age_wave:cPES_z", roi %in% tg) %>%
    transmute(model=tag, region=roi, beta=round(beta,5), se=round(se,5),
              conf.low=round(conf.low,5), conf.high=round(conf.high,5),
              p=signif(p,3), p_fdr_hemi=round(p_fdr_hemi,4), p_fdr_global=round(p_fdr_global,4))
}
unadj <- grab(file.path(WSB,"03_cPES_BRAIN/LMR_EXTRA_THICKNESS_cPES_FINAL/01_MAIN_MODEL/02_MAIN_all_terms_FDR_LH_RH.csv"),"UNADJUSTED")
adj   <- grab(file.path(WSB,"03b_cPES_BRAIN_ADJ_CBCL/LMR_EXTRA_THICKNESS_cPES_FINAL/01_MAIN_MODEL/02_MAIN_all_terms_FDR_LH_RH.csv"),"ADJ_baselineCBCL")
roi_tab <- bind_rows(unadj, adj)
cat("\n=== Per-ROI age×cPES: unadjusted vs adjusted ===\n"); print(as.data.frame(roi_tab), row.names=FALSE)
write_csv(roi_tab, file.path(OUT,"23_24_perROI_ageXcPES_unadj_vs_adj.csv"))
cat("\nSaved review-analysis CSVs to _REVIEW_ANALYSES\n")
