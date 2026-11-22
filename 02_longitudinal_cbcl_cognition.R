################################################################################
# 02_longitudinal_cbcl_cognition.R
# Cumulative psychosocial exposome score (cPES) study - BHRC.
# Linear mixed models of CBCL total and the four cognitive outcomes across waves (per-term estimates).
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
# cPES (no any_substance_use) -> Longitudinal validation
#   PART A: Cognition (Executive Function: wm_fl, ic_fl, tp_fl, ef_fl)
#   PART B: CBCL total (psychopathology) across waves W0-W2
# Clean, reproducible re-implementation of 2_LMR_FINAL_Cognition_cPES_EXTRA.R
#   - Data sources fixed to Windows paths:
#       * Cognition (Exec) rds
#       * Covariates + CBCL rebuilt from BHRC .sav (user-approved)
#       * cPES = POOLED out-of-fold score WITHOUT any_substance_use
#   - Model formulas kept IDENTICAL to the original script.
#   - Fixes original execution-order bugs (PATH_OUT / all_models_tidy).
#   - Audit notes (faithful to original code, NOT changed):
#       * cognition models use raw `age` (not centered)  [audit M3]
#       * cognition FDR pooled across ALL fixed terms of the 4 models [audit M4]
#       * CBCL LRT uses RAW cl_tot and DROPS site        [audit M2]
#       * CBCL simple-slope pairwise uses Bonferroni      [audit M1]
################################################################################

suppressWarnings(suppressMessages({
  library(haven); library(dplyr); library(tidyr); library(readr)
  library(janitor); library(purrr)
  library(lme4); library(lmerTest); library(broom.mixed)
  library(performance); library(emmeans); library(caret); library(ggplot2)
}))

set.seed(241094)

PATH_BHRC <- "<BHRC_INPUT_DIR>/Grassi-Oliveira_167BHRC_2023_12_09.sav"
PATH_cPES  <- "${cPES_OUT}/01_EXWAS_EN_cPES/POOLED/07_cpes_oof_predictions.csv"
PATH_COG  <- "<COGNITION_DIR>/Exec"
PATH_OUT  <- "${cPES_OUT}/02_LONGITUDINAL_COGNITION"
if (!dir.exists(PATH_OUT)) dir.create(PATH_OUT, recursive = TRUE)

summ <- c()  # summary log lines
addln <- function(...) { line <- paste0(...); summ <<- c(summ, line); cat(line, "\n") }

## ---- Shared inputs -----------------------------------------------------------
addln("=== Loading inputs ===")
BHRC_LONG <- read_sav(PATH_BHRC, encoding = "latin1") %>% clean_names() %>%
  filter(redcap_event_name %in% c("wave0_arm_1","wave1_arm_1","wave2_arm_1")) %>%
  transmute(ident = as.character(ident), redcap_event_name,
            cl_tot = as.numeric(cl_tot), age = as.numeric(age),
            gender = gender, site = site_original)

cPES <- read_csv(PATH_cPES, show_col_types = FALSE) %>%
  transmute(ident = as.character(ident), cPES_score_OOF_z = as.numeric(cPES_score_OOF_z))
addln("cPES n idents = ", nrow(cPES), " (WITHOUT any_substance_use)")

Execf <- readRDS(PATH_COG) %>% clean_names()
exec_df <- Execf %>%
  transmute(ident = as.character(ident), redcap_event_name,
            wm_fl, ic_fl, tp_fl, ef_fl)

################################################################################
# PART A — COGNITION
################################################################################
addln("\n=== PART A: Cognition ===")
exec_final <- BHRC_LONG %>%
  select(ident, redcap_event_name, age, site, gender) %>%
  inner_join(cPES, by = "ident") %>%
  inner_join(exec_df, by = c("ident","redcap_event_name")) %>%
  mutate(wave = case_when(redcap_event_name=="wave0_arm_1"~0,
                          redcap_event_name=="wave1_arm_1"~1,
                          redcap_event_name=="wave2_arm_1"~2, TRUE~NA_real_),
         gender = as.factor(gender), site = as.factor(site))

dat_exec <- exec_final %>%
  mutate(time = as.numeric(wave),
         age_c = age - mean(age, na.rm = TRUE),
         env_c = cPES_score_OOF_z - mean(cPES_score_OOF_z, na.rm = TRUE))

outcome_vars <- c("wm_fl","ic_fl","tp_fl","ef_fl")
base_formula <- "~ time * env_c + gender + site + age + (1 + time | ident)"   # raw age (faithful)

list_of_models <- setNames(outcome_vars, outcome_vars) %>%
  map(~ lmer(as.formula(paste(.x, base_formula)), data = dat_exec, REML = FALSE))

# FIX (original bug): define tidy table BEFORE applying FDR
all_models_tidy <- list_of_models %>%
  imap_dfr(~ broom.mixed::tidy(.x, effect = "fixed") %>% mutate(Outcome = .y))

all_models_tidy_corrected <- all_models_tidy %>%
  mutate(p.adjusted.fdr = p.adjust(p.value, method = "fdr"))   # FDR across ALL terms of ALL models (faithful)

write_csv(all_models_tidy_corrected, file.path(PATH_OUT, "all_models_tidy_corrected.csv"))

# Headline: cPES (env_c) main effect on each cognitive outcome
cog_cpes <- all_models_tidy_corrected %>% filter(term == "env_c") %>%
  select(Outcome, estimate, std.error, statistic, p.value, p.adjusted.fdr)
addln("cPES (env_c) main effect on cognition:")
for (i in seq_len(nrow(cog_cpes))) addln(sprintf("  %-6s beta=%+.4f  SE=%.4f  t=%.2f  p=%.4f  FDR=%.4f",
      cog_cpes$Outcome[i], cog_cpes$estimate[i], cog_cpes$std.error[i],
      cog_cpes$statistic[i], cog_cpes$p.value[i], cog_cpes$p.adjusted.fdr[i]))

################################################################################
# PART B — CBCL longitudinal
################################################################################
addln("\n=== PART B: CBCL longitudinal ===")
cbcl_df <- BHRC_LONG %>% select(ident, redcap_event_name, cl_tot, age, site, gender)
cbcl_m <- cPES %>% inner_join(cbcl_df, by = "ident")

pre_clt <- caret::preProcess(cbcl_m["cl_tot"], method = c("YeoJohnson","center","scale"))
cbcl_m$cl_tot_tr <- predict(pre_clt, cbcl_m["cl_tot"])[[1]]

dat <- cbcl_m %>%
  mutate(wave = case_when(redcap_event_name=="wave0_arm_1"~0,
                          redcap_event_name=="wave1_arm_1"~1,
                          redcap_event_name=="wave2_arm_1"~2, TRUE~NA_real_),
         gender = as.factor(gender), site = as.factor(site),
         time = as.numeric(wave),
         env_c = cPES_score_OOF_z - mean(cPES_score_OOF_z, na.rm = TRUE),
         age_c = age - mean(age, na.rm = TRUE))

mod <- lmer(cl_tot_tr ~ time * env_c + age_c + gender + site + (1 + time | ident),
            data = dat, REML = FALSE)

# Fixed effects (lmerTest t/df/p) + Wald CI
coefs <- as.data.frame(summary(mod)$coefficients)
coefs$term <- rownames(coefs)
ci <- as.data.frame(confint(mod, method = "Wald"))
coefs <- coefs %>% left_join(tibble(term = rownames(ci), CI_low = ci[[1]], CI_high = ci[[2]]), by = "term")
write_csv(coefs, file.path(PATH_OUT, "fixed_effects_lmm_cbcl.csv"))

gv <- function(t, col) coefs[coefs$term==t, col]
addln(sprintf("time      beta=%+.4f SE=%.4f df=%.0f t=%.2f p=%.3g",
      gv("time","Estimate"), gv("time","Std. Error"), gv("time","df"), gv("time","t value"), gv("time","Pr(>|t|)")))
addln(sprintf("cPES(env_c) beta=%+.4f SE=%.4f df=%.0f t=%.2f p=%.3g",
      gv("env_c","Estimate"), gv("env_c","Std. Error"), gv("env_c","df"), gv("env_c","t value"), gv("env_c","Pr(>|t|)")))
addln(sprintf("time:env_c beta=%+.4f SE=%.4f df=%.0f t=%.2f p=%.3g",
      gv("time:env_c","Estimate"), gv("time:env_c","Std. Error"), gv("time:env_c","df"), gv("time:env_c","t value"), gv("time:env_c","Pr(>|t|)")))

# Nakagawa R2 + ICC
r2n <- performance::r2_nakagawa(mod)
ic  <- performance::icc(mod)
addln(sprintf("Marginal R2 = %.4f | Conditional R2 = %.4f", r2n$R2_marginal, r2n$R2_conditional))
addln(sprintf("ICC adjusted = %.4f | ICC unadjusted = %.4f", ic$ICC_adjusted, ic$ICC_unadjusted))

# LRT (faithful: RAW cl_tot, NO site)
mod0 <- lmer(cl_tot ~ time + age_c + gender + (1 + time | ident), data = dat, REML = FALSE)
mod1 <- lmer(cl_tot ~ time * env_c + age_c + gender + (1 + time | ident), data = dat, REML = FALSE)
lrt <- anova(mod0, mod1)
addln(sprintf("LRT mod0 vs mod1: Chisq=%.2f df=%d p=%.3g",
      lrt$Chisq[2], lrt$Df[2], lrt$`Pr(>Chisq)`[2]))

# Simple slopes of time at cPES -1SD/0/+1SD (Bonferroni pairwise, faithful)
SDz <- sd(dat$cPES_score_OOF_z, na.rm = TRUE)
slopes <- emtrends(mod, specs = ~ env_c, var = "time", at = list(env_c = c(-SDz, 0, SDz)))
slopes_df <- as.data.frame(slopes)
write_csv(slopes_df, file.path(PATH_OUT, "cbcl_simple_slopes_time.csv"))
pw <- as.data.frame(pairs(slopes, adjust = "bonferroni"))
write_csv(pw, file.path(PATH_OUT, "cbcl_simple_slopes_pairwise_bonf.csv"))
addln("Simple slopes of time by cPES level:")
for (i in seq_len(nrow(slopes_df))) addln(sprintf("  env_c=%+.3f  slope=%+.4f  SE=%.4f",
      slopes_df$env_c[i], slopes_df$time.trend[i], slopes_df$SE[i]))

# Trajectory figure (main)
env_vals <- c(-SDz, 0, SDz)
pred_grid <- expand.grid(time = c(0,1,2), env_c = env_vals, age_c = 0,
                         gender = factor(1, levels=levels(dat$gender)),
                         site = factor(levels(dat$site)[1], levels=levels(dat$site)))
pred_grid$pred <- predict(mod, newdata = pred_grid, re.form = NA)
pred_grid$env_label <- factor(ifelse(pred_grid$env_c<0,"Low (-1SD)",ifelse(pred_grid$env_c>0,"High (+1SD)","Mean")),
                              levels=c("Low (-1SD)","Mean","High (+1SD)"))
gg <- ggplot(pred_grid, aes(time, pred, color=env_label)) +
  geom_line(linewidth=1.2) + geom_point(size=3) +
  scale_color_manual(values=c("Low (-1SD)"="darkgreen","Mean"="orange","High (+1SD)"="darkred"), name="cPES") +
  scale_x_continuous(breaks=c(0,1,2), labels=c("W0","W1","W2")) +
  labs(x="Wave", y="Predicted CBCL (transformed)", title="CBCL trajectories by cPES (no any_substance_use)") +
  theme_bw(14) + theme(legend.position="bottom")
ggsave(file.path(PATH_OUT, "Figure1_Trajectories_cPES.png"), gg, width=8, height=6, dpi=300)

writeLines(summ, file.path(PATH_OUT, "SUMMARY_longitudinal.txt"))
addln("\n=== DONE. Outputs in ", PATH_OUT, " ===")
