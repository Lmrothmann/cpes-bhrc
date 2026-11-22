################################################################################
# 06_allmetrics_sweep.R
# Cumulative psychosocial exposome score (cPES) study - BHRC.
# age x cPES across area/volume/curvature/folding/subcortical metrics (thickness-specificity check).
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
# Reviewer #14/#15: age×cPES across ALL FreeSurfer metrics, ComBat-harmonized.
# Metrics: thickness, area, volume, meancurv, foldind, curvind (cortical DK, LH+RH)
#          + aseg_volume (subcortical). Male-reference age×cPES + sex interaction.
################################################################################
suppressWarnings(suppressMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
  library(lme4); library(lmerTest); library(broom.mixed); library(longCombat)
}))
set.seed(241094)
TAB <- "<PROJECT_ROOT>/Output/tables"
WSB <- "${cPES_OUT}"
OUT <- file.path(WSB, "_REVIEW_ANALYSES"); if(!dir.exists(OUT)) dir.create(OUT, recursive=TRUE)

cPES <- read_csv(file.path(WSB,"01_EXWAS_EN_cPES/POOLED/07_cpes_oof_predictions.csv"), show_col_types=FALSE) %>%
  transmute(ident=as.character(ident), cPES_z=as.numeric(scale(cPES_score_OOF_z)))
BH <- read_sav <- haven::read_sav(file.path("<BHRC_INPUT_DIR>/Grassi-Oliveira_167BHRC_2023_12_09.sav"), encoding="latin1") %>% janitor::clean_names()
COV <- BH %>% transmute(ident=as.character(ident),
                        wave=dplyr::case_when(redcap_event_name=="wave0_arm_1"~0L,redcap_event_name=="wave1_arm_1"~1L,redcap_event_name=="wave2_arm_1"~2L,TRUE~NA_integer_),
                        age_wave=as.numeric(age),
                        gender=factor(as.integer(gender),levels=c(1,2),labels=c("Male","Female")),
                        site=factor(site_original)) %>% filter(!is.na(wave))

parse_id <- function(df){ df %>% rename(fs_id=1) %>%
  mutate(ident=as.character(as.integer(str_extract(fs_id,"(?<=sub-)\\d+"))),
         wave=dplyr::case_when(str_extract(fs_id,"(?<=_)ses-\\d+")=="ses-1"~0L,
                               str_extract(fs_id,"(?<=_)ses-\\d+")=="ses-2"~1L,
                               str_extract(fs_id,"(?<=_)ses-\\d+")=="ses-3"~2L,TRUE~NA_integer_)) }

analyze <- function(tab){
  raw <- suppressMessages(read_tsv(file.path(TAB, paste0(tab,".tsv")), show_col_types=FALSE))
  d <- parse_id(raw) %>% inner_join(cPES,by="ident") %>% inner_join(COV,by=c("ident","wave")) %>%
    filter(!is.na(age_wave),!is.na(cPES_z),!is.na(gender),!is.na(site))
  feats <- setdiff(names(raw)[-1], names(COV))
  feats <- feats[sapply(d[feats], function(x) is.numeric(x) && var(x,na.rm=TRUE)>1e-8)]
  d <- d %>% filter(if_all(all_of(feats), ~ !is.na(.x)))
  if(nrow(d)<50 || nlevels(droplevels(d$site))<2) return(NULL)
  cb <- tryCatch(longCombat(idvar="ident",timevar="wave",batchvar="site",features=feats,
                            formula="age_wave*cPES_z*gender",ranef="(1|ident)",data=as.data.frame(d),verbose=FALSE),
                 error=function(e) NULL)
  if(is.null(cb)) return(NULL)
  H <- cb$data_combat; hc <- grep("\\.combat$",names(H),value=TRUE)
  dd <- bind_cols(d %>% select(ident,age_wave,cPES_z,gender), H[hc])
  res <- lapply(hc, function(cc){
    dd$value <- dd[[cc]]
    m <- tryCatch(lmer(value~age_wave*cPES_z*gender+(1|ident),dd,REML=FALSE,control=lmerControl(calc.derivs=FALSE)),error=function(e)NULL)
    if(is.null(m)) return(NULL)
    t <- broom.mixed::tidy(m,effects="fixed")
    male <- t[t$term=="age_wave:cPES_z",]; sxd <- t[t$term=="age_wave:cPES_z:genderFemale",]
    data.frame(table=tab, roi=sub("\\.combat$","",cc),
               male_b=male$estimate, male_p=male$p.value,
               sexdiff_b=ifelse(nrow(sxd),sxd$estimate,NA), sexdiff_p=ifelse(nrow(sxd),sxd$p.value,NA))
  }) %>% bind_rows()
  if(nrow(res)==0) return(NULL)
  res$male_fdr <- p.adjust(res$male_p,"BH"); res
}

TABLES <- c(outer(c("lh_aparc_","rh_aparc_"), c("thickness","area","volume","meancurv","foldind","curvind"), paste0), "aseg_volume")
all <- bind_rows(lapply(TABLES, function(t){ cat("...",t,"\n"); analyze(t) }))
write_csv(all, file.path(OUT,"14_15_ALLmetrics_postComBat_ageXcPES.csv"))

cat("\n=== Survivors (male age×cPES, per-table FDR<0.05) by metric ===\n")
print(all %>% group_by(table) %>% summarise(n_roi=n(), n_sig_FDR05=sum(male_fdr<0.05,na.rm=TRUE),
        n_sexmod=sum(sexdiff_p<0.05,na.rm=TRUE)) %>% as.data.frame(), row.names=FALSE)

cat("\n=== The 3 thickness target regions ACROSS metrics (male age×cPES) ===\n")
core <- c("lateralorbitofrontal","inferiorparietal","supramarginal")
print(all %>% filter(str_detect(roi, paste(core,collapse="|")), str_starts(table,"rh_")) %>%
        transmute(table=str_remove(table,"rh_aparc_"), roi=str_remove(str_remove(roi,"rh_"),"_.*"),
                  male_b=round(male_b,5), male_fdr=round(male_fdr,3), sexdiff_p=round(sexdiff_p,3)) %>% as.data.frame(), row.names=FALSE)
