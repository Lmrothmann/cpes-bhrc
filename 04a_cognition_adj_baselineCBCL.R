################################################################################
# 04a_cognition_adj_baselineCBCL.R
# Cumulative psychosocial exposome score (cPES) study - BHRC.
# Cognitive models re-fitted with baseline CBCL total as a covariate.
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
# Reviewer #16/#17: cognition models ADJUSTED for baseline CBCL
# Does cPES predict cognition beyond the baseline psychopathology it was built on?
################################################################################
suppressWarnings(suppressMessages({
  library(haven); library(dplyr); library(readr); library(janitor); library(purrr)
  library(lme4); library(lmerTest); library(broom.mixed)
}))
set.seed(241094)

PATH_BHRC <- "<BHRC_INPUT_DIR>/Grassi-Oliveira_167BHRC_2023_12_09.sav"
PATH_cPES  <- "${cPES_OUT}/01_EXWAS_EN_cPES/POOLED/07_cpes_oof_predictions.csv"
PATH_COG  <- "<COGNITION_DIR>/Exec"
OUT       <- "${cPES_OUT}/02b_COGNITION_ADJ_CBCL"
if (!dir.exists(OUT)) dir.create(OUT, recursive = TRUE)
summ <- c(); addln <- function(...) { l<-paste0(...); summ<<-c(summ,l); cat(l,"\n") }

B <- read_sav(PATH_BHRC, encoding="latin1") %>% clean_names()
covs <- B %>% filter(redcap_event_name %in% c("wave0_arm_1","wave1_arm_1","wave2_arm_1")) %>%
  transmute(ident=as.character(ident), redcap_event_name, age=as.numeric(age),
            gender=gender, site=site_original)
# baseline CBCL (W0) as person-level covariate
cbcl0 <- B %>% filter(redcap_event_name=="wave0_arm_1") %>%
  transmute(ident=as.character(ident), cbcl0=as.numeric(cl_tot)) %>% filter(!is.na(cbcl0))
cPES <- read_csv(PATH_cPES, show_col_types=FALSE) %>%
  transmute(ident=as.character(ident), cPES_score_OOF_z=as.numeric(cPES_score_OOF_z))
Ex <- readRDS(PATH_COG) %>% clean_names() %>%
  transmute(ident=as.character(ident), redcap_event_name, wm_fl, ic_fl, tp_fl, ef_fl)

dat <- covs %>% inner_join(cPES, by="ident") %>%
  inner_join(Ex, by=c("ident","redcap_event_name")) %>%
  inner_join(cbcl0, by="ident") %>%
  mutate(time = case_when(redcap_event_name=="wave0_arm_1"~0,
                          redcap_event_name=="wave1_arm_1"~1,
                          redcap_event_name=="wave2_arm_1"~2, TRUE~NA_real_),
         gender=as.factor(gender), site=as.factor(site),
         env_c = cPES_score_OOF_z - mean(cPES_score_OOF_z, na.rm=TRUE),
         cbcl0_c = cbcl0 - mean(cbcl0, na.rm=TRUE))

outcomes <- c("wm_fl","ic_fl","tp_fl","ef_fl")
# UNADJUSTED (reference) and ADJUSTED (+ cbcl0_c) — env_c main effect
fit_set <- function(adj){
  fml <- if (adj) "~ time * env_c + cbcl0_c + gender + site + age + (1 + time | ident)"
         else     "~ time * env_c + gender + site + age + (1 + time | ident)"
  setNames(outcomes, outcomes) %>%
    map(~ lmer(as.formula(paste(.x, fml)), data=dat, REML=FALSE)) %>%
    imap_dfr(~ broom.mixed::tidy(.x, effect="fixed") %>% mutate(Outcome=.y)) %>%
    mutate(FDR = p.adjust(p.value, "fdr"))
}
unadj <- fit_set(FALSE); adj <- fit_set(TRUE)

addln("=== cPES (env_c) main effect on cognition: UNADJUSTED vs ADJUSTED for baseline CBCL ===")
for (o in outcomes) {
  u <- unadj %>% filter(Outcome==o, term=="env_c"); a <- adj %>% filter(Outcome==o, term=="env_c")
  addln(sprintf("  %-6s UNADJ beta=%+.4f p=%.4f FDR=%.4f  |  ADJ beta=%+.4f p=%.4f FDR=%.4f",
    o, u$estimate, u$p.value, u$FDR, a$estimate, a$p.value, a$FDR))
}
# baseline-CBCL effect itself (for context)
addln("\n=== baseline-CBCL (cbcl0_c) effect in adjusted models ===")
for (o in outcomes) { c0 <- adj %>% filter(Outcome==o, term=="cbcl0_c")
  addln(sprintf("  %-6s cbcl0 beta=%+.5f p=%.4g", o, c0$estimate, c0$p.value)) }

write_csv(adj, file.path(OUT, "cognition_ADJ_baselineCBCL_all_terms.csv"))
write_csv(unadj, file.path(OUT, "cognition_UNADJ_all_terms.csv"))
writeLines(summ, file.path(OUT, "SUMMARY_cognition_adj.txt"))
addln("\nDONE ->", OUT)
