################################################################################
# 09_bootstrap_calibration.R
# Cumulative psychosocial exposome score (cPES) study - BHRC.
# Bootstrap 95% CIs and calibration for out-of-fold performance.
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
# Reviewer #19/#20: bootstrap 95% CIs for OOF metrics (R2/RMSE/rho) +
# calibration slope/intercept, pooled and sex-stratified.
################################################################################
suppressWarnings(suppressMessages({library(haven); library(dplyr); library(readr); library(janitor); library(caret)}))
set.seed(241094)

PATH_BHRC <- "<BHRC_INPUT_DIR>/Grassi-Oliveira_167BHRC_2023_12_09.sav"
PATH_EXPO <- "<PROJECT_ROOT>/Output/BHRCS_final_Exposome.sav"
WS  <- "${cPES_OUT}/01_EXWAS_EN_cPES"
OUT <- "${cPES_OUT}/_REVIEW_ANALYSES"
if(!dir.exists(OUT)) dir.create(OUT, recursive=TRUE)

# Rebuild analytic sample + outcome_tr exactly as script 1
BHRC <- read_sav(PATH_BHRC, encoding="latin1") %>% clean_names()
w0 <- BHRC %>% filter(redcap_event_name=="wave0_arm_1") %>% transmute(ident=as.character(ident), cl_tot=as.numeric(cl_tot), age=as.numeric(age), gender, site_original)
Expo <- read_sav(PATH_EXPO) %>% clean_names() %>% select(-any_of(c("cl_tot","age","gender","site_original")))
Expo$ident <- as.character(Expo$ident)
dat <- Expo %>% inner_join(w0, by="ident") %>%
  filter(!is.na(cl_tot), !is.na(age), !is.na(gender), !is.na(site_original))
preY <- caret::preProcess(dat["cl_tot"], method=c("YeoJohnson","center","scale"))
dat$outcome_tr <- predict(preY, dat["cl_tot"])[[1]]
dat <- dat %>% select(ident, gender, cl_tot, outcome_tr)

metrics <- function(pred, obs_tr, obs_raw){
  ok <- !is.na(pred) & !is.na(obs_tr)
  R2 <- suppressWarnings(cor(pred[ok], obs_tr[ok]))^2
  RMSE <- sqrt(mean((pred[ok]-obs_tr[ok])^2))
  rho <- suppressWarnings(cor(pred[ok], obs_raw[ok], method="spearman"))
  cal <- lm(obs_tr ~ pred); ci <- suppressWarnings(confint(cal))
  c(R2=R2, RMSE=RMSE, rho=rho,
    cal_int=coef(cal)[1], cal_slope=coef(cal)[2],
    cal_int_lo=ci[1,1], cal_int_hi=ci[1,2], cal_slope_lo=ci[2,1], cal_slope_hi=ci[2,2])
}
bootci <- function(df, predcol, B=2000){
  n <- nrow(df); est <- metrics(df[[predcol]], df$outcome_tr, df$cl_tot)
  bs <- matrix(NA, B, 3)
  for(b in 1:B){ i <- sample.int(n, n, replace=TRUE); s <- df[i,]
    ok <- !is.na(s[[predcol]]); p <- s[[predcol]][ok]; ot <- s$outcome_tr[ok]; or <- s$cl_tot[ok]
    bs[b,1] <- suppressWarnings(cor(p,ot))^2; bs[b,2] <- sqrt(mean((p-ot)^2)); bs[b,3] <- suppressWarnings(cor(p,or,method="spearman")) }
  q <- function(x) quantile(x, c(.025,.975), na.rm=TRUE)
  list(est=est, R2ci=q(bs[,1]), RMSEci=q(bs[,2]), rhoci=q(bs[,3]))
}

rows <- list()
for(grp in c("POOLED","MALE_1","FEMALE_2")){
  cpes <- read_csv(file.path(WS,grp,"07_cpes_oof_predictions.csv"), show_col_types=FALSE) %>% transmute(ident=as.character(ident), cPES_pred_OOF)
  en  <- read_csv(file.path(WS,grp,"02_enet_oof_predictions.csv"), show_col_types=FALSE) %>% transmute(ident=as.character(ident), EN_pred_OOF)
  d <- dat %>% inner_join(cpes, by="ident") %>% inner_join(en, by="ident")
  for(mod in c("cPES_pred_OOF","EN_pred_OOF")){
    r <- bootci(d, mod)
    rows[[paste(grp,mod)]] <- data.frame(
      stratum=grp, model=ifelse(mod=="cPES_pred_OOF","cPES","EN"), n=nrow(d),
      R2=round(r$est["R2"],3), R2_lo=round(r$R2ci[1],3), R2_hi=round(r$R2ci[2],3),
      RMSE=round(r$est["RMSE"],3), RMSE_lo=round(r$RMSEci[1],3), RMSE_hi=round(r$RMSEci[2],3),
      rho=round(r$est["rho"],3), rho_lo=round(r$rhoci[1],3), rho_hi=round(r$rhoci[2],3),
      cal_intercept=round(r$est["cal_int.(Intercept)"],3), cal_slope=round(r$est["cal_slope.pred"],3),
      cal_slope_lo=round(r$est["cal_slope_lo"],3), cal_slope_hi=round(r$est["cal_slope_hi"],3))
  }
}
res <- do.call(rbind, rows); rownames(res) <- NULL
write_csv(res, file.path(OUT,"19_20_bootstrap_calibration.csv"))
cat("=== Bootstrap CIs + calibration (2000 reps) ===\n"); print(res, row.names=FALSE)
cat("\nSaved: 19_20_bootstrap_calibration.csv\n")
