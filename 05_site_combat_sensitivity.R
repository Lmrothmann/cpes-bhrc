################################################################################
# 05_site_combat_sensitivity.R
# Cumulative psychosocial exposome score (cPES) study - BHRC.
# Site as fixed vs random effect and longitudinal ComBat, male reference. (Supersedes the buggy alphabetical-gender ComBat script.)
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
# CORRECTED: age×cPES (MALE reference, matching original script 3's reported term)
# under site-FIXED, site-RANDOM, and ComBat. Male ref = factor(levels Male,Female).
################################################################################
suppressWarnings(suppressMessages({library(dplyr); library(readr); library(lme4); library(lmerTest); library(broom.mixed); library(longCombat)}))
set.seed(241094)
WSB <- "${cPES_OUT}"
QCD <- file.path(WSB, "03_cPES_BRAIN/LMR_EXTRA_THICKNESS_cPES_FINAL/00_QC")
OUT <- file.path(WSB, "_REVIEW_ANALYSES")
TERM <- "age_wave:cPES_z"   # = age×cPES at MALE reference

axe <- function(m){ if(is.null(m)) return(c(NA,NA)); t<-broom.mixed::tidy(m,effects="fixed"); r<-t[t$term==TERM,]; if(nrow(r)==0) return(c(NA,NA)); c(r$estimate, r$p.value) }

hemi <- function(tag){
  wm <- read_csv(file.path(QCD, paste0("00_", tag, "_wide_merged.csv")), show_col_types=FALSE)
  rois <- grep("_thickness$", names(wm), value=TRUE); rois <- rois[!grepl("MeanThickness|meanthickness", rois)]
  d <- wm %>% mutate(ident=as.character(ident),
                     gender=factor(gender, levels=c("Male","Female")),  # MALE reference (as in script 3)
                     site=factor(site_original)) %>%
    filter(!is.na(age_wave),!is.na(cPES_z),!is.na(gender),!is.na(site)) %>%
    filter(if_all(all_of(rois), ~ !is.na(.x)))
  # ComBat harmonized features (preserve age*cPES*gender)
  cb <- longCombat(idvar="ident", timevar="wave", batchvar="site", features=rois,
                   formula="age_wave*cPES_z*gender", ranef="(1|ident)", data=as.data.frame(d), verbose=FALSE)
  harm <- cb$data_combat; hcols <- grep("\\.combat$", names(harm), value=TRUE)
  dC <- bind_cols(d %>% select(ident,wave,age_wave,cPES_z,gender), harm[hcols])

  out <- lapply(rois, function(cc){
    d$value <- as.numeric(d[[cc]])
    mf <- tryCatch(lmer(value~age_wave*cPES_z*gender+site+(1|ident), d, REML=FALSE, control=lmerControl(calc.derivs=FALSE)), error=function(e) NULL)
    mr <- tryCatch(suppressWarnings(lmer(value~age_wave*cPES_z*gender+(1|ident)+(1|site), d, REML=FALSE, control=lmerControl(calc.derivs=FALSE))), error=function(e) NULL)
    dC$value <- dC[[paste0(cc,".combat")]]
    mc <- tryCatch(lmer(value~age_wave*cPES_z*gender+(1|ident), dC, REML=FALSE, control=lmerControl(calc.derivs=FALSE)), error=function(e) NULL)
    data.frame(roi=cc, fx_b=axe(mf)[1], fx_p=axe(mf)[2], rd_b=axe(mr)[1], rd_p=axe(mr)[2], cb_b=axe(mc)[1], cb_p=axe(mc)[2])
  }) %>% bind_rows()
  out$fx_fdr <- p.adjust(out$fx_p,"BH"); out$rd_fdr <- p.adjust(out$rd_p,"BH"); out$cb_fdr <- p.adjust(out$cb_p,"BH")
  out$hemi <- tag; out
}
all <- bind_rows(hemi("lh_aparc_thickness"), hemi("rh_aparc_thickness"))
write_csv(all, file.path(OUT,"CORRECTED_site_combat_maleRef.csv"))
tg <- c("rh_lateralorbitofrontal_thickness","rh_inferiorparietal_thickness","rh_supramarginal_thickness")
cat("=== age×cPES (MALE ref) in 3 target ROIs: site-FIXED | site-RANDOM | ComBat ===\n")
print(all %>% filter(roi %in% tg) %>% transmute(roi=sub("_thickness","",roi),
  fixed_b=round(fx_b,5), fixed_fdr=round(fx_fdr,3),
  random_b=round(rd_b,5), random_fdr=round(rd_fdr,3),
  combat_b=round(cb_b,5), combat_fdr=round(cb_fdr,3)), row.names=FALSE)
cat(sprintf("\nROIs surviving per-hemi FDR<0.05:  site-fixed=%d | site-random=%d | ComBat=%d\n",
  sum(all$fx_fdr<0.05,na.rm=TRUE), sum(all$rd_fdr<0.05,na.rm=TRUE), sum(all$cb_fdr<0.05,na.rm=TRUE)))
