################################################################################
# 01_build_cPES.R
# Cumulative psychosocial exposome score (cPES) study - BHRC.
# Nested-CV elastic net; out-of-fold cPES; full-sample coefficients (eTable S2); ExWAS; selection stability.
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
# Exposome Approach - BHRCS 
# - HARD DROP before NA handling
# - Complete only on outcome + covariates
# - ExWAS (LM | LMM) with global median imputation (no leakage)
# - EN/cPES with per-fold imputation (train median) — no leakage
# - cPES/EN OOF saved per individual (row_id + ident)
# - Plots prontos p/ paper
# Author: Leonardo Melo Rothmann - 2025
# (contact via corresponding author)
################################################################################

set.seed(241094)
options(stringsAsFactors = FALSE)

# --------------------------- 0) Caminhos e pacotes ----------------------------
PATH_OUT  <- "<PROJECT_ROOT>/SENSITIVITY_SUBSTANCE/BASELINE_RERUN"
PATH_BHRC <- "<BHRC_INPUT_DIR>/Grassi-Oliveira_167BHRC_2023_12_09.sav"
PATH_EXPO <- "<PROJECT_ROOT>/Output/BHRCS_final_Exposome.sav"
if (!dir.exists(PATH_OUT)) dir.create(PATH_OUT, recursive = TRUE)

# How to handle 'site_original': "covariate" (default) | "random" (ExWAS uses lmer) | "none"
SITE_MODE       <- "covariate"
RUN_POOLED      <- TRUE
RUN_BY_SEX      <- TRUE
DROP_CBCL_ITEMS <- TRUE
SAVE_QQPLOT     <- TRUE
FAST_MODE       <- FALSE        # TRUE = runs faster (fewer grids/repeats)

libs <- c("tidyverse","haven","janitor","broom","readr","caret","glmnet","ggpubr","patchwork",
          "lme4","broom.mixed")
for (p in libs) { if (!requireNamespace(p, quietly=TRUE)) message(">> Install: ", p) }
invisible(lapply(intersect(libs, rownames(installed.packages())), library, character.only=TRUE))

# Tema/cores p/ paper
COL_RISK <- "#D55E00"; COL_PROT <- "#0072B2"; COL_NEUTRO <- "grey40"
theme_pub <- theme_bw(12) +
  theme(panel.grid.major = element_line(color="grey90", linewidth=.25),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face="bold", size=13),
        legend.position = "right")

# salvar PNG + PDF (com fallback)
ggsave_both <- function(path_png, plot_obj, w=8, h=6, dpi=300){
  ggplot2::ggsave(path_png, plot_obj, width=w, height=h, dpi=dpi)
  # try PDF; if Cairo is not available, use "pdf"
  ok <- TRUE
  pdf_path <- sub("\\.png$", ".pdf", path_png)
  tryCatch({
    ggplot2::ggsave(pdf_path, plot_obj, width=w, height=h, device=grDevices::cairo_pdf)
  }, error = function(e){ ok <<- FALSE })
  if (!ok) ggplot2::ggsave(pdf_path, plot_obj, width=w, height=h, device="pdf")
}

# per-group title helper
group_title <- function(label){
  if (grepl("^POOLED$", label)) "All sample"
  else if (grepl("^MALE", label)) "Male"
  else if (grepl("^FEMALE", label)) "Female"
  else label
}

# --------------------------- 1) Ler e montar baseline -------------------------
BHRC_all <- haven::read_sav(PATH_BHRC, encoding="latin1") %>% janitor::clean_names()
stopifnot(all(c("ident","redcap_event_name") %in% names(BHRC_all)))
keep_cols <- intersect(c("ident","redcap_event_name","cl_tot","age","gender","site_original"),
                       names(BHRC_all))

BHRC_keep <- BHRC_all %>%
  dplyr::select(dplyr::all_of(keep_cols)) %>%
  dplyr::mutate(wave = dplyr::case_when(
    redcap_event_name == "wave0_arm_1" ~ 0L,
    redcap_event_name == "wave1_arm_1" ~ 1L,
    redcap_event_name == "wave2_arm_1" ~ 2L,
    TRUE ~ NA_integer_
  )) %>% dplyr::filter(!is.na(wave))

w0 <- BHRC_keep %>% dplyr::filter(wave==0L) %>% dplyr::select(ident, cl_tot, age)

Expo <- haven::read_sav(PATH_EXPO) %>% janitor::clean_names()
stopifnot("ident" %in% names(Expo))
Expo_clean <- Expo %>% dplyr::select(-dplyr::any_of(c("cl_tot","age")))
dat_base0 <- Expo_clean %>% dplyr::inner_join(w0, by="ident")

# candidate covariates
COVARS_CAND <- intersect(c("age","gender","site_original"), names(dat_base0))
if (SITE_MODE == "none") COVARS_CAND <- setdiff(COVARS_CAND, "site_original")

# tmp (pre-imputation) — keep ident!
ALL_EXPO0 <- setdiff(names(dat_base0), c("ident","cl_tot", COVARS_CAND))
tmp <- dat_base0 %>%
  dplyr::transmute(ident = ident,
                   outcome = as.numeric(cl_tot),
                   !!!rlang::syms(COVARS_CAND),
                   !!!rlang::syms(ALL_EXPO0))

# obvious cleanups
for (v in c("reg36b1","reg36b2","reg36b3","reg36b4","reg36b5","reg36b6",
            "trauma_exp","traumac_exp","traumap_exp","frelation")) {
  if (v %in% names(tmp)) tmp[[v]] <- NULL
}

# --------------------------- 2) Missing (completo em Y+covars) ----------------
dat_complete <- tmp %>%
  dplyr::filter(!is.na(outcome)) %>%
  dplyr::filter(if_all(dplyr::all_of(COVARS_CAND), ~ !is.na(.x)))

# --------------------------- 3) Transformar Y / derivar features --------------
preY <- caret::preProcess(dat_complete["outcome"], method=c("YeoJohnson","center","scale"))
dat_complete$outcome_tr <- predict(preY, dat_complete["outcome"])[[1]]

# parental SUD: 1 if ANY {alcdep, alcabus, drugdep, drugabus}==1; NA if all NA
v_sud <- intersect(c("imini_cur_alcdep","imini_cur_alcabus","imini_cur_drugdep","imini_cur_drugabus"),
                   names(dat_complete))
if (length(v_sud) > 0) {
  mat1 <- sapply(v_sud, function(nm) as.integer(dat_complete[[nm]]==1))
  all_na <- apply(sapply(v_sud, function(nm) is.na(dat_complete[[nm]])), 1, all)
  any_1  <- apply(mat1, 1, function(r) any(r==1, na.rm=TRUE))
  dat_complete$imini_cur_sud <- ifelse(all_na, NA_integer_, ifelse(any_1, 1L, 0L))
}

dat <- dat_complete %>%
  dplyr::mutate(gender = as.factor(gender),
                site_original = as.factor(site_original))
dat$row_id <- seq_len(nrow(dat))

# final covariates
BASE_COVARS <- intersect(c("age","gender","site_original"), names(dat))
if (SITE_MODE == "none") BASE_COVARS <- setdiff(BASE_COVARS, "site_original")

# candidate exposures
ADMIN_VARS <- c("respondent","ident")
ALL_EXPO <- setdiff(names(dat), c("ident","outcome","outcome_tr", BASE_COVARS, "row_id"))
ALL_EXPO <- setdiff(ALL_EXPO, intersect(ADMIN_VARS, names(dat)))

# Remove CBCL-derived? (off by default)
CBCL_PATTERNS <- c("^cl\\d+$", "^cbcl")
if (DROP_CBCL_ITEMS) {
  cbcl_drop <- unique(unlist(lapply(CBCL_PATTERNS, function(p) grep(p, names(dat), ignore.case=TRUE, value=TRUE))))
  cbcl_drop <- setdiff(cbcl_drop, "cl_tot")
  ALL_EXPO <- setdiff(ALL_EXPO, cbcl_drop)
  readr::write_csv(tibble::tibble(dropped_cbcl = cbcl_drop), file.path(PATH_OUT,"_dropped_cbcl_from_X.csv"))
}

# numeric only
is_num <- vapply(dat[, ALL_EXPO, drop=FALSE], function(x) is.numeric(x) && !is.factor(x), logical(1))
ALL_EXPO <- ALL_EXPO[is_num]

# --------------------------- 3b) HARD DROP (adjustments/default) -------------------
DROP_HARD <- c(
  "cat_mother_ageatbirth","maecrianca","maetem_hs","abepstrat",
  "prenatal_alcohol","prenatal_smooking",
  "neonatal_icu","prematurely",
  "cl2","cl39","cl55","cl99","cl105","bad_cia",
  "p_bulliedhim","p_bulliedhim_fr","bullyingbother",
  "father_known",
  "pphyab","pphyneg","psexab","pemoab",
  "imini_cur_any","imini_cur_any_mood","imini_cur_any_anx","imini_cur_any_subs",
  "divorce","fat20a_y","mrelation","reg36","reg38"
)
DROP_HARD <- c(DROP_HARD, intersect(c("imini_cur_drugdep","imini_cur_drugabus",
                                      "imini_cur_alcdep","imini_cur_alcabus"), names(dat)))

ALL_EXPO <- setdiff(ALL_EXPO, DROP_HARD)

# keep only the necessary columns in dat, no duplicates
cols_keep <- unique(c("ident","outcome","outcome_tr", BASE_COVARS, "row_id", ALL_EXPO, "imini_cur_sud"))
dat <- dat[, intersect(cols_keep, names(dat)), drop = FALSE]
dat <- dat[, !duplicated(names(dat)), drop = FALSE]

# rebuild ALL_EXPO from the final dat
ALL_EXPO <- setdiff(names(dat), c("ident","outcome","outcome_tr", BASE_COVARS, "row_id"))

# Zero-variance (drop only zeroVar)
if (length(ALL_EXPO) > 0) {
  nz <- caret::nearZeroVar(dat[, ALL_EXPO, drop=FALSE], saveMetrics=TRUE)
  zero_drop <- rownames(nz)[nz$zeroVar]
  if (length(zero_drop)) for (v in zero_drop) dat[[v]] <- NULL
}
ALL_EXPO <- intersect(ALL_EXPO, names(dat))

message("N = ", nrow(dat),
        " | p(exposures) = ", length(ALL_EXPO),
        " | covars = ", paste(BASE_COVARS, collapse=", "))

# --------------------------- 4) Utils/plots -----------------------------------
qqplot_p <- function(pvals, path_png){
  chis  <- qchisq(pmax(1e-300, 1 - pvals), df=1)
  lambda <- median(chis, na.rm=TRUE)/qchisq(0.5, df=1)
  if (SAVE_QQPLOT) {
    obs  <- -log10(sort(pvals))
    exp  <- -log10(ppoints(sum(!is.na(pvals))))
    png(path_png, w=800, h=800, res=120)
    plot(exp, obs, xlab="Expected -log10(p)", ylab="Observed -log10(p)", main="ExWAS QQ")
    abline(0,1)
    dev.off()
  }
  invisible(lambda)
}

# --------------------------- 5) EXWAS (LM/LMM) --------------------------------
run_exwas <- function(dat0, exposures, covars, outdir, grp_title = "All sample"){
  # global median imputation
  X_raw <- dat0[, exposures, drop=FALSE]
  meds  <- sapply(X_raw, function(x) median(x, na.rm=TRUE))
  for (nm in names(X_raw)) if (anyNA(X_raw[[nm]])) X_raw[[nm]][is.na(X_raw[[nm]])] <- meds[[nm]]
  preX <- caret::preProcess(X_raw, method=c("YeoJohnson","center","scale"))
  Xz   <- predict(preX, X_raw)
  use_lmm <- (SITE_MODE == "random") && ("site_original" %in% names(dat0))
  base_cols <- unique(c("outcome","outcome_tr", covars, if (use_lmm) "site_original"))
  datX <- cbind(dat0[, base_cols, drop=FALSE], Xz)
  
  fit_one <- function(v){
    if (use_lmm) {
      covs_no_site <- setdiff(covars, "site_original")
      f <- as.formula(paste("outcome_tr ~", paste(c(v, covs_no_site), collapse=" + "), "+ (1|site_original)"))
      fit <- try(lme4::lmer(f, data=datX, REML=FALSE), silent=TRUE)
      if (inherits(fit,"try-error")) return(data.frame(variable=v, beta=NA, se=NA, t=NA, p=NA))
      tt <- try(broom.mixed::tidy(fit), silent=TRUE)
      if (inherits(tt,"try-error")) return(data.frame(variable=v, beta=NA, se=NA, t=NA, p=NA))
      row <- tt[tt$effect=="fixed" & tt$term==v, c("estimate","std.error","statistic","p.value")]
      if (nrow(row)==0) return(data.frame(variable=v, beta=NA, se=NA, t=NA, p=NA))
      data.frame(variable=v, beta=row$estimate, se=row$std.error, t=row$statistic, p=row$p.value)
    } else {
      f <- as.formula(paste("outcome_tr ~", paste(c(v, covars), collapse=" + ")))
      fit <- try(lm(f, data=datX), silent=TRUE)
      if (inherits(fit,"try-error")) return(data.frame(variable=v, beta=NA, se=NA, t=NA, p=NA))
      tt <- broom::tidy(fit)
      row <- tt[tt$term==v, c("estimate","std.error","statistic","p.value")]
      if (nrow(row)==0) return(data.frame(variable=v, beta=NA, se=NA, t=NA, p=NA))
      data.frame(variable=v, beta=row$estimate, se=row$std.error, t=row$statistic, p=row$p.value)
    }
  }
  
  exwas <- dplyr::bind_rows(lapply(exposures, fit_one))
  exwas$pFDR <- p.adjust(exwas$p, "BH")
  exwas$direction <- ifelse(is.na(exwas$beta), NA, ifelse(exwas$beta>=0,"risk","protective"))
  readr::write_csv(exwas, file.path(outdir,"01_exwas_results.csv"))
  
  # Forest with ALL (dynamic height)
  ex_all <- exwas %>% arrange(desc(abs(beta))) %>%
    mutate(variable = factor(variable, levels = rev(variable)))
  n_for <- nrow(ex_all); h_for <- max(6, min(0.18*n_for + 2, 40))
  p_ex <- ggplot(ex_all, aes(x=variable, y=beta, color=direction)) +
    geom_point() +
    geom_errorbar(aes(ymin=beta-1.96*se, ymax=beta+1.96*se), width=.2) +
    coord_flip() + theme_pub +
    labs(x="Exposures (z)",
         y=paste0("β (", ifelse(use_lmm,"LMM","LM"), ", outcome_tr)"),
         title=paste0("ExWAS — ", grp_title)) +
    scale_color_manual(values=c("risk"=COL_RISK,"protective"=COL_PROT))
  ggsave_both(file.path(outdir,"01_exwas_forest.png"), p_ex, w=8, h=h_for)
  
  # Volcano (per-group title)
  exwas$neglog10 <- -log10(pmax(exwas$pFDR, 1e-300))
  exwas$sig <- ifelse(exwas$pFDR<=0.05,"FDR≤0.05","NS")
  p_vol <- ggplot(exwas, aes(x=beta, y=neglog10, color=direction, shape=sig))+
    geom_hline(yintercept=-log10(0.05), linetype=2, color="grey50")+
    geom_point(alpha=.9, size=2)+ theme_pub +
    labs(x=paste0("β (", ifelse(use_lmm,"LMM","LM"), ", outcome_tr)"),
         y="-log10(FDR)", title=paste0("ExWAS Volcano — ", grp_title))+
    scale_color_manual(values=c("risk"=COL_RISK,"protective"=COL_PROT))
  ggsave_both(file.path(outdir,"01_exwas_volcano.png"), p_vol, w=8.5, h=6)
  
  lambda <- qqplot_p(exwas$p, file.path(outdir,"01_exwas_qqplot.png"))
  writeLines(paste("Lambda_GC ~", round(lambda,3)), con=file.path(outdir,"01_lambda_gc.txt"))
  invisible(exwas)
}

# --------------------------- 6) Folds -----------------------------------------
make_folds <- function(datY, k=10){
  bins <- cut(datY$outcome_tr, breaks=quantile(datY$outcome_tr, seq(0,1,0.2), na.rm=TRUE), include.lowest=TRUE)
  if (!("site_original" %in% names(datY))) {
    strata <- factor(bins)
  } else {
    strata <- interaction(datY$site_original, bins, drop=TRUE)
  }
  caret::createFolds(factor(strata), k=k, list=TRUE)
}

# --------------------------- 7) EN OOF ----------------------------------------
run_enet_oof <- function(dat0, exposures, outdir, grp_title = "All sample"){
  y <- dat0$outcome_tr
  folds <- make_folds(dat0, k=10)
  pred_outer <- rep(NA_real_, nrow(dat0))
  sel_list <- list(); oof_cor_raw <- rep(NA_real_, length(folds))
  
  tl <- if (FAST_MODE) 10 else 40
  ctrl_in <- caret::trainControl(method="repeatedcv", number=10, repeats=if (FAST_MODE) 2 else 5)
  
  for (k in seq_along(folds)) {
    te <- folds[[k]]; tr <- setdiff(seq_len(nrow(dat0)), te)
    vars_tr <- exposures[vapply(dat0[tr, exposures, drop=FALSE],
                                function(x) var(x, na.rm=TRUE) > 0, logical(1))]
    if (length(vars_tr)==0) { pred_outer[te] <- mean(y[tr]); next }
    
    X_tr_raw <- dat0[tr, vars_tr, drop=FALSE]
    med_tr   <- sapply(X_tr_raw, function(x) median(x, na.rm=TRUE))
    Xtr_imp  <- as.data.frame(Map(function(x,m) ifelse(is.na(x), m, x), X_tr_raw, med_tr))
    X_te_raw <- dat0[te, vars_tr, drop=FALSE]
    Xte_imp  <- as.data.frame(Map(function(x,m) ifelse(is.na(x), m, x), X_te_raw, med_tr))
    
    preX_tr <- caret::preProcess(Xtr_imp, method=c("YeoJohnson","center","scale"))
    Xtr <- predict(preX_tr, Xtr_imp) %>% as.matrix()
    Xte <- predict(preX_tr, Xte_imp) %>% as.matrix()
    ytr <- y[tr]
    
    en_k <- caret::train(x=Xtr, y=ytr, method="glmnet", trControl=ctrl_in, tuneLength=tl)
    pred_outer[te] <- as.numeric(predict(en_k, newdata=Xte))
    
    cf_k <- as.matrix(coef(en_k$finalModel, en_k$bestTune$lambda))[-1,,drop=FALSE]
    nz <- which(as.numeric(cf_k)!=0)
    if (length(nz)>0) sel_list[[length(sel_list)+1]] <- data.frame(variable=rownames(cf_k)[nz], fold=k)
    
    oof_cor_raw[k] <- suppressWarnings(cor(pred_outer[te], dat0$outcome[te], method="spearman", use="pairwise.complete.obs"))
  }
  
  RMSE_en <- caret::RMSE(pred_outer, y)
  R2_en   <- caret::R2(pred_outer, y)
  rho_en  <- suppressWarnings(cor(pred_outer, dat0$outcome, method="spearman", use="pairwise.complete.obs"))
  perf_en <- tibble::tibble(Model="Elastic Net (nestedCV, no leakage)",
                            RMSE=RMSE_en, R2=R2_en, Spearman_with_rawY=rho_en,
                            MeanFoldSpearman_rawY = mean(oof_cor_raw, na.rm=TRUE))
  readr::write_csv(perf_en, file.path(outdir,"02_enet_performance_oof.csv"))
  
  # OOF per individual
  readr::write_csv(tibble::tibble(ident = dat0$ident, row_id = dat0$row_id,
                                  EN_pred_OOF = pred_outer),
                   file.path(outdir,"02_enet_oof_predictions.csv"))
  
  # calibration with per-group title
  df_cal <- tibble::tibble(obs=dat0$outcome_tr, pred=pred_outer)
  p_cal <- ggplot(df_cal, aes(x=obs, y=pred))+
    geom_point(alpha=.3, size=1, color="#555555")+
    geom_smooth(method="lm", se=FALSE, linewidth=.8, color=COL_RISK)+
    geom_abline(slope=1, intercept=0, linetype=2, color="grey50")+
    annotate("text", x=min(df_cal$obs,na.rm=TRUE), y=max(df_cal$pred,na.rm=TRUE), hjust=0, vjust=1,
             label=paste0("R²=", sprintf("%.3f", R2_en), " | RMSE=", sprintf("%.3f", RMSE_en),
                          " | ρ=", sprintf("%.3f", rho_en)))+
    theme_pub + labs(x="Observed (outcome_tr)", y="Predicted (EN, OOF)",
                     title=paste0("Calibration — Elastic Net (OOF) — ", grp_title))
  ggsave_both(file.path(outdir,"02_enet_oof_calibration.png"), p_cal, w=7.5, h=6)
  
  # Selection frequency — ALL
  if (length(sel_list)) {
    sel_df <- dplyr::bind_rows(sel_list) %>% dplyr::count(variable, name="n_folds_selected") %>%
      dplyr::mutate(freq = n_folds_selected/length(folds)) %>% dplyr::arrange(dplyr::desc(freq))
    readr::write_csv(sel_df, file.path(outdir,"02_enet_selection_freq.csv"))
    sel_all <- sel_df %>% mutate(variable = factor(variable, levels=rev(variable)))
    n_sel <- nrow(sel_all); h_sel <- max(6, min(0.18*n_sel + 2, 40))
    p_sel <- ggplot(sel_all, aes(x=variable, y=freq))+
      geom_col() + coord_flip()+ theme_pub +
      labs(x="Exposure", y="Selection frequency (10-fold)",
           title=paste0("Elastic Net — Stability — ", grp_title))
    ggsave_both(file.path(outdir,"02_enet_selection_freq.png"), p_sel, w=7.5, h=h_sel)
  }
  
  # Final fit (interpretation; NOT OOF) — saves “06_…”
  X_all_raw <- dat0[, exposures, drop=FALSE]
  med_all   <- sapply(X_all_raw, function(x) median(x, na.rm=TRUE))
  X_all_imp <- as.data.frame(Map(function(x,m) ifelse(is.na(x), m, x), X_all_raw, med_all))
  preX_all  <- caret::preProcess(X_all_imp, method=c("YeoJohnson","center","scale"))
  X_all     <- predict(preX_all, X_all_imp) %>% as.matrix()
  ctrl_final <- caret::trainControl(method="repeatedcv", number=10, repeats=if (FAST_MODE) 2 else 5)
  enet_final <- caret::train(x=X_all, y=y, method="glmnet", trControl=ctrl_final, tuneLength=tl)
  cf_final <- as.matrix(coef(enet_final$finalModel, enet_final$bestTune$lambda))[-1,,drop=FALSE]
  coef_enet <- tibble::tibble(variable=rownames(cf_final), beta=as.numeric(cf_final)) %>% dplyr::filter(beta!=0)
  
  # saves compatibility "02_" and new "06_"
  readr::write_csv(coef_enet, file.path(outdir,"02_enet_nonzero_coefs_FULLDATA_INTERPRET_ONLY.csv"))
  readr::write_csv(coef_enet, file.path(outdir,"06_enet_nonzero_coefs.csv"))
  
  if (nrow(coef_enet) > 0) {
    coef_plot <- coef_enet %>%
      mutate(effect = if_else(beta >= 0, "risk", "protective")) %>%
      arrange(desc(abs(beta))) %>%
      mutate(variable = factor(variable, levels = rev(variable)))
    n_coef <- nrow(coef_plot); h_coef <- max(6, min(0.20*n_coef + 1.5, 40))
    p_en <- ggplot(coef_plot, aes(x = variable, y = beta, color = effect)) +
      geom_point() +
      geom_segment(aes(xend = variable, y = 0, yend = beta), linewidth = .6) +
      coord_flip() + theme_bw(12) +
      labs(x = "Exposures (z)", y = "β (EN, outcome_tr)",
           title = paste0("Elastic Net — ", grp_title)) +
      scale_color_manual(values = c("risk"=COL_RISK,"protective"=COL_PROT))
    ggsave(file.path(outdir,"06_enet_coef.png"), p_en, width=8, height=h_coef, dpi=300)
    ggsave(file.path(outdir,"06_enet_coef.pdf"), p_en, width=8, height=h_coef, device="pdf")
  }
  
  invisible(list(pred_oof=pred_outer, coef_enet=coef_enet))
}

# --------------------------- 8) cPES OOF ---------------------------------------
run_cpes_oof <- function(dat0, exposures, covars, outdir){
  folds <- make_folds(dat0, k=10)
  pred_cpes <- rep(NA_real_, nrow(dat0))
  cpes_score_oof <- rep(NA_real_, nrow(dat0))
  tl <- if (FAST_MODE) 10 else 40
  
  for (k in seq_along(folds)) {
    te <- folds[[k]]; tr <- setdiff(seq_len(nrow(dat0)), te)
    vars_tr <- exposures[vapply(dat0[tr, exposures, drop=FALSE],
                                function(x) var(x, na.rm=TRUE) > 0, logical(1))]
    if (length(vars_tr)==0) { 
      pred_cpes[te] <- mean(dat0$outcome_tr[tr]); 
      cpes_score_oof[te] <- 0
      next 
    }
    
    X_tr_raw <- dat0[tr, vars_tr, drop=FALSE]
    med_tr   <- sapply(X_tr_raw, function(x) median(x, na.rm=TRUE))
    Xtr_imp  <- as.data.frame(Map(function(x,m) ifelse(is.na(x), m, x), X_tr_raw, med_tr))
    X_te_raw <- dat0[te, vars_tr, drop=FALSE]
    Xte_imp  <- as.data.frame(Map(function(x,m) ifelse(is.na(x), m, x), X_te_raw, med_tr))
    
    preX_tr <- caret::preProcess(Xtr_imp, method=c("YeoJohnson","center","scale"))
    Xtr <- predict(preX_tr, Xtr_imp) %>% as.matrix()
    Xte <- predict(preX_tr, Xte_imp) %>% as.matrix()
    ytr <- dat0$outcome_tr[tr]
    
    ctrl_in <- caret::trainControl(method="repeatedcv", number=10, repeats=if (FAST_MODE) 2 else 5)
    en_k <- caret::train(x=Xtr, y=ytr, method="glmnet", trControl=ctrl_in, tuneLength=tl)
    cf_k <- as.matrix(coef(en_k$finalModel, en_k$bestTune$lambda))[-1,,drop=FALSE]
    w <- as.numeric(cf_k); names(w) <- rownames(cf_k); w <- w[w!=0]
    if (length(w)==0) { 
      pred_cpes[te] <- mean(ytr); 
      cpes_score_oof[te] <- 0
      next 
    }
    
    cPES_tr_raw <- as.numeric(Xtr[, names(w), drop=FALSE] %*% w)
    sgn <- ifelse(suppressWarnings(cor(cPES_tr_raw, dat0$outcome[tr], method="spearman", use="pairwise.complete.obs")) < 0, -1, 1)
    w <- w * sgn
    cPES_te <- as.numeric(Xte[, names(w), drop=FALSE] %*% w)
    
    cov_use <- if (SITE_MODE=="covariate") covars else setdiff(covars, "site_original")
    df_tr <- data.frame(y=ytr, cPES=cPES_tr_raw*sgn, dat0[tr, cov_use, drop=FALSE])
    df_te <- data.frame(cPES=cPES_te,      dat0[te, cov_use, drop=FALSE])
    m_k <- lm(as.formula(paste("y ~ cPES +", paste(cov_use, collapse=" + "))), data=df_tr)
    pred_cpes[te] <- as.numeric(predict(m_k, newdata=df_te))
    cpes_score_oof[te] <- cPES_te
  }
  
  RMSE_cpes <- caret::RMSE(pred_cpes, dat0$outcome_tr)
  R2_cpes   <- caret::R2(pred_cpes, dat0$outcome_tr)
  rho_cpes  <- suppressWarnings(cor(pred_cpes, dat0$outcome, method="spearman", use="pairwise.complete.obs"))
  perf_cpes <- tibble::tibble(Model="cPES_EN + covariates (proper OOF, no leakage)",
                             RMSE=RMSE_cpes, R2=R2_cpes, Spearman_with_rawY=rho_cpes)
  readr::write_csv(perf_cpes, file.path(outdir,"07_perf_cPES_EN_oof.csv"))
  
  cpes_z <- as.numeric(scale(cpes_score_oof))
  readr::write_csv(tibble::tibble(ident = dat0$ident, row_id = dat0$row_id,
                                  cPES_score_OOF = cpes_score_oof,
                                  cPES_score_OOF_z = cpes_z,
                                  cPES_pred_OOF = pred_cpes),
                   file.path(outdir,"07_cpes_oof_predictions.csv"))
  
  df_cal <- tibble::tibble(obs=dat0$outcome_tr, pred=pred_cpes)
  p_cal <- ggplot(df_cal, aes(x=obs, y=pred))+
    geom_point(alpha=.3, size=1, color="#555555")+
    geom_smooth(method="lm", se=FALSE, linewidth=.8, color=COL_PROT)+
    geom_abline(slope=1, intercept=0, linetype=2, color="grey50")+
    annotate("text", x=min(df_cal$obs,na.rm=TRUE), y=max(df_cal$pred,na.rm=TRUE), hjust=0, vjust=1,
             label=paste0("R²=", sprintf("%.3f", R2_cpes), " | RMSE=", sprintf("%.3f", RMSE_cpes),
                          " | ρ=", sprintf("%.3f", rho_cpes)))+
    theme_pub + labs(x="Observed (outcome_tr)", y="Predicted (cPES, OOF)", title="Calibration — cPES (OOF)")
  ggsave_both(file.path(outdir,"07_cpes_oof_calibration.png"), p_cal, w=7.5, h=6)
  
  invisible(pred_cpes)
}

# --------------------------- 9) Site sensitivity (optional) -----------------
exwas_site_sensitivity <- function(datZ, exposures, covars_base, outdir){
  out_with <- file.path(outdir,"_tmp_with_site"); out_no <- file.path(outdir,"_tmp_no_site")
  dir.create(out_with, showWarnings=FALSE, recursive=TRUE)
  dir.create(out_no,   showWarnings=FALSE, recursive=TRUE)
  
  exwas_with <- run_exwas(datZ, exposures, covars=covars_base, outdir=out_with, grp_title="All sample")
  exwas_no   <- run_exwas(datZ, exposures, covars=setdiff(covars_base,"site_original"), outdir=out_no, grp_title="All sample")
  
  unlink(out_with, recursive=TRUE); unlink(out_no, recursive=TRUE)
  cmp <- exwas_with %>% dplyr::select(variable, beta_with=beta, pFDR_with=pFDR) %>%
    dplyr::left_join(exwas_no %>% dplyr::select(variable, beta_no=beta, pFDR_no=pFDR), by="variable") %>%
    dplyr::mutate(delta_beta = beta_with - beta_no, flip_dir = sign(beta_with)!=sign(beta_no))
  readr::write_csv(cmp, file.path(outdir,"01_exwas_compare_site_vs_nosite.csv"))
}

# --------------------------- 10) Runner (POOLED/SEX) --------------------------
run_pipeline <- function(dat_full, exposures, covars, label){
  outdir <- file.path(PATH_OUT, label); if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
  grp <- group_title(label)
  
  message("\n========== PIPELINE: ", label, " ==========")
  message("N = ", nrow(dat_full), " | covars = ", paste(covars, collapse=", "),
          " | p = ", length(exposures), " | SITE_MODE = ", SITE_MODE)
  
  # 10.1 ExWAS
  exwas  <- run_exwas(dat_full, exposures, covars, outdir, grp_title = grp)
  
  # 10.2 EN nested OOF + 06_*
  en_obj <- run_enet_oof(dat_full, exposures, outdir, grp_title = grp)
  
  # 10.3 cPES OOF
  covars_use <- if (SITE_MODE=="covariate") covars else setdiff(covars, "site_original")
  run_cpes_oof(dat_full, exposures, covars_use, outdir)
  
  # 10.4 Master
  tb <- exwas %>% dplyr::select(variable, exwas_beta=beta, exwas_se=se, exwas_p=p, exwas_pFDR=pFDR)
  if (!is.null(en_obj$coef_enet) && nrow(en_obj$coef_enet)>0) {
    tb <- tb %>% dplyr::left_join(en_obj$coef_enet %>% dplyr::select(variable, en_beta=beta), by="variable")
  }
  tb <- tb %>% dplyr::arrange(dplyr::desc(abs(dplyr::coalesce(en_beta, exwas_beta))))
  readr::write_csv(tb, file.path(outdir,"09_master_table_exposome.csv"))
  
  if ("site_original" %in% names(dat_full) && SITE_MODE=="covariate") {
    exwas_site_sensitivity(dat_full, exposures, covars, outdir)
  }
  
  sink(file.path(outdir,"sessionInfo.txt")); print(sessionInfo()); sink()
  message(">> Done:", normalizePath(outdir))
}

# --------------------------- 11) Rodadas --------------------------------------
which_dup <- which(duplicated(names(dat))); if (length(which_dup)) { print(names(dat)[which_dup]) }

dat_pooled <- dat %>% dplyr::mutate(gender=factor(gender))
COVARS_POOLED <- intersect(c("age","gender", if (SITE_MODE=="covariate") "site_original"), names(dat_pooled))
if (SITE_MODE=="none") COVARS_POOLED <- setdiff(COVARS_POOLED, "site_original")

if (RUN_POOLED) run_pipeline(dat_pooled, ALL_EXPO, COVARS_POOLED, "POOLED")

if (RUN_BY_SEX && "gender" %in% names(dat_pooled)) {
  lvl <- levels(dat_pooled$gender)
  if (length(lvl) >= 2) {
    dat_m <- dat_pooled %>% dplyr::filter(gender == lvl[1]) %>% droplevels()
    COVARS_M <- intersect(c("age","site_original"), names(dat_m))
    if (SITE_MODE=="none") COVARS_M <- setdiff(COVARS_M, "site_original")
    if (nrow(dat_m) > 50) run_pipeline(dat_m, ALL_EXPO, COVARS_M, paste0("MALE_", lvl[1])) else message("Small male N:", nrow(dat_m))
    
    dat_f <- dat_pooled %>% dplyr::filter(gender == lvl[2]) %>% droplevels()
    COVARS_F <- intersect(c("age","site_original"), names(dat_f))
    if (SITE_MODE=="none") COVARS_F <- setdiff(COVARS_F, "site_original")
    if (nrow(dat_f) > 50) run_pipeline(dat_f, ALL_EXPO, COVARS_F, paste0("FEMALE_", lvl[2])) else message("Small female N:", nrow(dat_f))
  } else {
    message(">> 'gender' with <2 levels; stratification not possible.")
  }
}

message(">> Overall end. Outputs in:", normalizePath(PATH_OUT))
################################################################################
################################################################################
################################################################################