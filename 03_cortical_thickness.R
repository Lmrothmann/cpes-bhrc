################################################################################
# 03_cortical_thickness.R
# Cumulative psychosocial exposome score (cPES) study - BHRC.
# age x cPES x sex mixed models per Desikan-Killiany parcel (per-term estimates). Male is the reference level.
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
# FS THICKNESS + cPES — FINAL LONGITUDINAL PIPELINE (ONE SCRIPT, LH -> RH -> MERGE)
#
# PRIMARY MODEL
#   value ~ age_wave * cPES_z * gender + site + (1 | ident)
#
# ROBUSTNESS MODEL
#   value ~ age_within * cPES_z * gender + age_mean_c + site + (1 | ident)
#
# WORKFLOW
#   STEP 1: Run LH
#   STEP 2: Run RH
#   STEP 3: Merge LH + RH
#   STEP 4: Compute FDR (hemi-specific and global)
#   STEP 5: Generate unified tables and figures
#
# NOTE
#   - Plots are generated for ALL shortlisted ROIs, not only top hits
#   - Hemisphere-level FDR is used as the primary threshold for figures
################################################################################

set.seed(241094)
options(stringsAsFactors = FALSE)

suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(janitor)
  library(haven)
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(ggplot2)
  library(ggeffects)
  library(purrr)
  library(forcats)
  library(tibble)
})

HAS_GGSEG <- requireNamespace("ggseg", quietly = TRUE)

# ==============================================================================
# 0) PATHS / SETTINGS
# ==============================================================================
setwd("<PROJECT_ROOT>")

PATH_FS_TABLES <- "<PROJECT_ROOT>/Output/tables"
PATH_cPES       <- "${cPES_OUT}/01_EXWAS_EN_cPES/POOLED/07_cpes_oof_predictions.csv"
PATH_BHRC      <- "<BHRC_INPUT_DIR>/Grassi-Oliveira_167BHRC_2023_12_09.sav"

OUT_BASE <- "${cPES_OUT}/03_cPES_BRAIN/LMR_EXTRA_THICKNESS_cPES_FINAL"
OUT_QC   <- file.path(OUT_BASE, "00_QC")
OUT_MAIN <- file.path(OUT_BASE, "01_MAIN_MODEL")
OUT_ROB  <- file.path(OUT_BASE, "02_ROBUSTNESS_WITHIN_BETWEEN")
OUT_TAB  <- file.path(OUT_BASE, "03_TABLES_MANUSCRIPT")
OUT_FIGS <- file.path(OUT_BASE, "04_FIGURES")
OUT_MAPS <- file.path(OUT_BASE, "05_BRAINMAP_DATA")
OUT_TMP  <- file.path(OUT_BASE, "99_HEMI_INTERMEDIATE")

dir.create(OUT_QC, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_MAIN, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_ROB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_FIGS, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_MAPS, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TMP, recursive = TRUE, showWarnings = FALSE)

USE_SITE <- TRUE
FDR_MAIN <- 0.05
FDR_FALLBACK <- 0.10

AGE_MIN_OK <- 5
AGE_MAX_OK <- 27

DPI_FIG <- 600
TOP_PER_TERM_FOREST <- 20

TERM_MAIN_cPES      <- "cPES_z"
TERM_MAIN_AXE      <- "age_wave:cPES_z"
TERM_MAIN_AXE_SEX  <- "age_wave:cPES_z:genderFemale"

TERM_ROB_cPES       <- "cPES_z"
TERM_ROB_AXE       <- "age_within:cPES_z"
TERM_ROB_AXE_SEX   <- "age_within:cPES_z:genderFemale"

KEY_TERMS_MAIN <- c(TERM_MAIN_cPES, TERM_MAIN_AXE, TERM_MAIN_AXE_SEX)
KEY_TERMS_ROB  <- c(TERM_ROB_cPES, TERM_ROB_AXE, TERM_ROB_AXE_SEX)

# visual constants
COL_LOW  <- "#2c7bb6"  # azul
COL_HIGH <- "#d7191c"  # vermelho
COL_MEAN <- "black"

# ==============================================================================
# 1) HELPERS
# ==============================================================================
msg <- function(...) message(sprintf(...))

pick_existing <- function(cands, nm){
  out <- cands[cands %in% nm][1]
  if (is.na(out)) return(NA_character_)
  out
}

safe_nlevels <- function(x){
  if (is.factor(x)) return(nlevels(x))
  length(unique(x))
}

parse_fs_id <- function(df, id_col = 1){
  df %>%
    rename(fs_id = !!id_col) %>%
    mutate(
      ident = as.character(as.integer(str_extract(fs_id, "(?<=sub-)\\d+"))),
      ses   = str_extract(fs_id, "(?<=_)ses-\\d+"),
      wave  = case_when(
        ses == "ses-1" ~ 0L,
        ses == "ses-2" ~ 1L,
        ses == "ses-3" ~ 2L,
        TRUE ~ NA_integer_
      )
    ) %>%
    select(-ses)
}

table_meta <- function(table_name){
  hemi <- ifelse(str_starts(table_name, "lh_"), "LH",
                 ifelse(str_starts(table_name, "rh_"), "RH", NA))
  atlas <- ifelse(str_detect(table_name, "aparc"), "aparc", "unknown")
  tibble(table = table_name, hemi = hemi, atlas = atlas)
}

term_label_main <- function(term){
  dplyr::case_when(
    term == TERM_MAIN_cPES     ~ "Main cPES",
    term == TERM_MAIN_AXE     ~ "Age × cPES",
    term == TERM_MAIN_AXE_SEX ~ "Age × cPES × Female",
    TRUE ~ term
  )
}

term_label_rob <- function(term){
  dplyr::case_when(
    term == TERM_ROB_cPES      ~ "Main cPES",
    term == TERM_ROB_AXE      ~ "Within-age × cPES",
    term == TERM_ROB_AXE_SEX  ~ "Within-age × cPES × Female",
    TRUE ~ term
  )
}

theme_paper <- function(base_size = 12){
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = rel(0.95)),
      axis.title = element_text(face = "plain"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.3)
    )
}


NON_CORTICAL_COLS <- c(
  "etiv", "e_tiv", "brain_seg_vol_not_vent", "brainsegvolnotvent",
  
  "lh_mean_thickness_thickness", "rh_mean_thickness_thickness",
  "lh_mean_thickness", "rh_mean_thickness",
  
  "lh_meanthickness_thickness", "rh_meanthickness_thickness",
  "lh_meanthickness", "rh_meanthickness",
  
  "mean_thickness", "meanthickness_thickness",
  "meanthickness", "mean"
)

is_cortical_roi <- function(x){
  xx <- tolower(x)
  
  is_hemi <- grepl("^(lh|rh)_", xx)
  is_thick <- grepl("_thickness$", xx)
  
  is_non_cortical <- xx %in% NON_CORTICAL_COLS |
    grepl("mean[_]?thickness", xx) |
    grepl("^lh_mean", xx) |
    grepl("^rh_mean", xx)
  
  is_hemi & is_thick & !is_non_cortical
}

resolve_roi_col <- function(wide_df, roi_pick){
  nm <- names(wide_df)
  if (roi_pick %in% nm) return(roi_pick)
  
  cand <- roi_pick
  cand2 <- cand %>% str_replace("^lh_", "") %>% str_replace("^rh_", "")
  cand3 <- cand2 %>%
    str_replace("_thickness_thickness$", "_thickness") %>%
    str_replace("_thickness$", "")
  
  try_exact <- unique(c(
    cand, cand2, cand3,
    paste0("lh_", cand2), paste0("rh_", cand2),
    paste0("lh_", cand3), paste0("rh_", cand3),
    paste0(cand2, "_thickness"),
    paste0(cand3, "_thickness")
  ))
  
  hit <- try_exact[try_exact %in% nm][1]
  if (!is.na(hit)) return(hit)
  
  core <- cand3 %>%
    str_replace("^lh_", "") %>%
    str_replace("^rh_", "") %>%
    str_replace("_thickness$", "")
  
  hits <- nm[str_detect(tolower(nm), fixed(tolower(core)))]
  hits2 <- hits[str_detect(tolower(hits), "thickness")]
  if (length(hits2) > 0) hits <- hits2
  
  if (length(hits) >= 1) {
    msg("ROI resolver: '%s' -> chose '%s'", roi_pick, hits[1])
    return(hits[1])
  }
  
  stop("Could not resolve ROI column for roi_pick='", roi_pick, "'.")
}

qc_wave_coverage <- function(df_wide_m, table_name){
  sub_wave_dist <- df_wide_m %>%
    distinct(ident, wave) %>%
    count(ident, name = "n_waves") %>%
    count(n_waves, name = "n_subjects") %>%
    arrange(n_waves) %>%
    mutate(table = table_name) %>%
    relocate(table)
  
  age_by_wave <- df_wide_m %>%
    group_by(wave) %>%
    summarise(
      table = table_name,
      n_rows = n(),
      n_sub = n_distinct(ident),
      age_min = min(age_wave, na.rm = TRUE),
      age_med = median(age_wave, na.rm = TRUE),
      age_max = max(age_wave, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(sub_wave_dist = sub_wave_dist, age_by_wave = age_by_wave)
}

pretty_roi_label <- function(roi){
  roi %>%
    str_remove("^lh_") %>%
    str_remove("^rh_") %>%
    str_remove("_thickness$") %>%
    str_replace_all("_", " ") %>%
    str_to_title()
}

get_hemi_from_roi <- function(roi){
  case_when(
    str_starts(roi, "lh_") ~ "left",
    str_starts(roi, "rh_") ~ "right",
    TRUE ~ NA_character_
  )
}

get_region_from_roi <- function(roi){
  roi %>%
    str_remove("^lh_") %>%
    str_remove("^rh_") %>%
    str_remove("_thickness$")
}

to_ggseg_region <- function(roi){
  core <- get_region_from_roi(roi) %>% tolower()
  
  dplyr::case_when(
    core == "bankssts" ~ "banksts",
    core == "caudalanteriorcingulate" ~ "caudal anterior cingulate",
    core == "caudalmiddlefrontal" ~ "caudal middle frontal",
    core == "cuneus" ~ "cuneus",
    core == "entorhinal" ~ "entorhinal",
    core == "frontalpole" ~ "frontal pole",
    core == "fusiform" ~ "fusiform",
    core == "inferiorparietal" ~ "inferior parietal",
    core == "inferiortemporal" ~ "inferior temporal",
    core == "insula" ~ "insula",
    core == "isthmuscingulate" ~ "isthmus cingulate",
    core == "lateraloccipital" ~ "lateral occipital",
    core == "lateralorbitofrontal" ~ "lateral orbitofrontal",
    core == "lingual" ~ "lingual",
    core == "medialorbitofrontal" ~ "medial orbitofrontal",
    core == "middletemporal" ~ "middle temporal",
    core == "parahippocampal" ~ "parahippocampal",
    core == "paracentral" ~ "paracentral",
    core == "parsopercularis" ~ "pars opercularis",
    core == "parsorbitalis" ~ "pars orbitalis",
    core == "parstriangularis" ~ "pars triangularis",
    core == "pericalcarine" ~ "pericalcarine",
    core == "postcentral" ~ "postcentral",
    core == "posteriorcingulate" ~ "posterior cingulate",
    core == "precentral" ~ "precentral",
    core == "precuneus" ~ "precuneus",
    core == "rostralanteriorcingulate" ~ "rostral anterior cingulate",
    core == "rostralmiddlefrontal" ~ "rostral middle frontal",
    core == "superiorfrontal" ~ "superior frontal",
    core == "superiorparietal" ~ "superior parietal",
    core == "superiortemporal" ~ "superior temporal",
    core == "supramarginal" ~ "supramarginal",
    core == "temporalpole" ~ "temporal pole",
    core == "transversetemporal" ~ "transverse temporal",
    core == "unknown" ~ "unknown",
    TRUE ~ core
  )
}

pick_reference_site <- function(x){
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

get_cpes_col_after_clean <- function(df){
  pick_existing(c("cpes_z", "cPES_z"), names(df))
}

# ==============================================================================
# 2) READ cPES + BHRC LONG
# ==============================================================================
msg("\n[1/10] Reading cPES ...")
cpes_raw <- read_csv(PATH_cPES, show_col_types = FALSE) %>% clean_names()
stopifnot("ident" %in% names(cpes_raw))

cpes_z_col <- pick_existing(c("cpes_score_oof_z", "cpes_score_oof_z_"), names(cpes_raw))
if (is.na(cpes_z_col)) stop("cPES z column not found in cPES file.")

cpes <- cpes_raw %>%
  transmute(
    ident = as.character(ident),
    cPES_score_OOF_z = as.numeric(.data[[cpes_z_col]])
  ) %>%
  filter(!is.na(cPES_score_OOF_z)) %>%
  distinct() %>%
  mutate(
    cPES_z = as.numeric(scale(cPES_score_OOF_z))
  )

msg("  cPES: n=%d unique idents", nrow(cpes))

msg("\n[2/10] Reading BHRC longitudinal covariates ...")
BHRC_all <- read_sav(PATH_BHRC, encoding = "latin1") %>% clean_names()
stopifnot(all(c("ident", "redcap_event_name") %in% names(BHRC_all)))

BHRC_long <- BHRC_all %>%
  transmute(
    ident = as.character(ident),
    redcap_event_name = as.character(redcap_event_name),
    wave = case_when(
      redcap_event_name == "wave0_arm_1" ~ 0L,
      redcap_event_name == "wave1_arm_1" ~ 1L,
      redcap_event_name == "wave2_arm_1" ~ 2L,
      TRUE ~ NA_integer_
    ),
    age_wave = suppressWarnings(as.numeric(age)),
    gender = factor(
      suppressWarnings(as.integer(gender)),
      levels = c(1, 2),
      labels = c("Male", "Female")
    ),
    site_original = factor(site_original)
  ) %>%
  filter(!is.na(wave))

msg("  BHRC_long rows: n=%d", nrow(BHRC_long))

# ==============================================================================
# 3) MODEL FUNCTIONS
# ==============================================================================
fit_main_roi <- function(dat_roi, use_site = TRUE){
  
  f <- as.formula(paste0(
    "value ~ age_wave * cPES_z * gender",
    if (use_site) " + site_original" else "",
    " + (1 | ident)"
  ))
  
  m <- tryCatch(
    lmer(f, data = dat_roi, REML = FALSE,
         control = lmerControl(calc.derivs = FALSE)),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  
  tidy(m, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      model = "MAIN_age_wave",
      model_random = "RI_only",
      singular = lme4::isSingular(m, tol = 1e-4),
      converged = is.null(m@optinfo$conv$lme4$messages)
    )
}

fit_robust_roi <- function(dat_roi, use_site = TRUE){
  
  f <- as.formula(paste0(
    "value ~ age_within * cPES_z * gender + age_mean_c",
    if (use_site) " + site_original" else "",
    " + (1 | ident)"
  ))
  
  m <- tryCatch(
    lmer(f, data = dat_roi, REML = FALSE,
         control = lmerControl(calc.derivs = FALSE)),
    error = function(e) NULL
  )
  if (is.null(m)) return(NULL)
  
  tidy(m, effects = "fixed", conf.int = TRUE) %>%
    mutate(
      model = "ROBUST_within_between",
      model_random = "RI_only",
      singular = lme4::isSingular(m, tol = 1e-4),
      converged = is.null(m@optinfo$conv$lme4$messages)
    )
}

# ==============================================================================
# 4) FUNCTION TO RUN ONE HEMISPHERE
# ==============================================================================
run_one_hemi <- function(file_tsv){
  
  table_name <- tools::file_path_sans_ext(basename(file_tsv))
  meta <- table_meta(table_name)
  msg("\n[3/10] Running hemisphere/table: %s", table_name)
  
  wide_raw <- read_tsv(file_tsv, show_col_types = FALSE)
  wide <- parse_fs_id(wide_raw, 1)
  
  wide_m <- wide %>%
    left_join(cpes, by = "ident") %>%
    left_join(
      BHRC_long %>% select(ident, wave, age_wave, gender, site_original),
      by = c("ident", "wave")
    ) %>%
    filter(!is.na(wave), !is.na(age_wave), !is.na(gender), !is.na(cPES_z))
  
  age_flag <- wide_m %>%
    summarise(
      table = table_name,
      n_rows = n(),
      n_sub = n_distinct(ident),
      age_min = min(age_wave, na.rm = TRUE),
      age_max = max(age_wave, na.rm = TRUE),
      n_age_lt_min = sum(age_wave < AGE_MIN_OK, na.rm = TRUE),
      n_age_gt_max = sum(age_wave > AGE_MAX_OK, na.rm = TRUE)
    )
  
  age_components <- wide_m %>%
    group_by(ident) %>%
    summarise(age_mean = mean(age_wave, na.rm = TRUE), .groups = "drop")
  
  grand_age_mean <- mean(age_components$age_mean, na.rm = TRUE)
  
  wide_m <- wide_m %>%
    left_join(age_components, by = "ident") %>%
    mutate(
      age_within = age_wave - age_mean,
      age_mean_c = age_mean - grand_age_mean
    )
  
  qc <- wide_m %>%
    summarise(
      table = table_name,
      hemi = meta$hemi[1],
      atlas = meta$atlas[1],
      n_rows = n(),
      n_sub = n_distinct(ident),
      n_wave0 = sum(wave == 0),
      n_wave1 = sum(wave == 1),
      n_wave2 = sum(wave == 2),
      n_ge2 = {
        tmp <- distinct(pick(everything()), ident, wave) %>% count(ident, name = "n_waves")
        sum(tmp$n_waves >= 2)
      },
      n_ge3 = {
        tmp <- distinct(pick(everything()), ident, wave) %>% count(ident, name = "n_waves")
        sum(tmp$n_waves >= 3)
      },
      site_levels = safe_nlevels(site_original),
      age_min = min(age_wave, na.rm = TRUE),
      age_max = max(age_wave, na.rm = TRUE)
    )
  
  qc2 <- qc_wave_coverage(wide_m, table_name)
  
  write_csv(wide_m, file.path(OUT_QC, paste0("00_", table_name, "_wide_merged.csv")))
  write_csv(qc, file.path(OUT_QC, paste0("00_", table_name, "_QC_overview.csv")))
  write_csv(age_flag, file.path(OUT_QC, paste0("00_", table_name, "_QC_age_flags.csv")))
  write_csv(qc2$sub_wave_dist, file.path(OUT_QC, paste0("00_", table_name, "_QC_sub_wave_dist.csv")))
  write_csv(qc2$age_by_wave, file.path(OUT_QC, paste0("00_", table_name, "_QC_age_by_wave.csv")))
  
  long_roi <- wide_m %>%
    pivot_longer(
      cols = -c(fs_id, ident, wave, cPES_score_OOF_z, cPES_z, age_wave, gender,
                site_original, age_mean, age_within, age_mean_c),
      names_to = "roi",
      values_to = "value"
    ) %>%
    mutate(value = suppressWarnings(as.numeric(value))) %>%
    filter(!is.na(value)) %>%
    mutate(table = table_name, hemi = meta$hemi[1], atlas = meta$atlas[1])
  
  use_site_now <- USE_SITE && (nlevels(long_roi$site_original) >= 2)
  rois <- sort(unique(long_roi$roi))
  
  sample_qc <- long_roi %>%
    group_by(table, hemi, atlas, roi) %>%
    summarise(
      n_obs = n(),
      n_sub = n_distinct(ident),
      n_wave0 = sum(wave == 0),
      n_wave1 = sum(wave == 1),
      n_wave2 = sum(wave == 2),
      n_ge2 = {
        tmp <- distinct(pick(everything()), ident, wave) %>% count(ident, name = "n_waves")
        sum(tmp$n_waves >= 2)
      },
      n_ge3 = {
        tmp <- distinct(pick(everything()), ident, wave) %>% count(ident, name = "n_waves")
        sum(tmp$n_waves >= 3)
      },
      value_mean = mean(value, na.rm = TRUE),
      value_sd = sd(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  res_main <- vector("list", length(rois))
  res_rob  <- vector("list", length(rois))
  
  for (i in seq_along(rois)) {
    r <- rois[i]
    dat_roi <- long_roi %>% filter(roi == r)
    
    m1 <- fit_main_roi(dat_roi, use_site = use_site_now)
    if (!is.null(m1)) {
      res_main[[i]] <- m1 %>%
        mutate(
          table = table_name,
          hemi = meta$hemi[1],
          atlas = meta$atlas[1],
          roi = r,
          n_obs = nrow(dat_roi),
          n_sub = n_distinct(dat_roi$ident),
          n_ge2 = {
            tmp <- dat_roi %>% distinct(ident, wave) %>% count(ident, name = "n_waves")
            sum(tmp$n_waves >= 2)
          },
          n_ge3 = {
            tmp <- dat_roi %>% distinct(ident, wave) %>% count(ident, name = "n_waves")
            sum(tmp$n_waves >= 3)
          }
        )
    }
    
    m2 <- fit_robust_roi(dat_roi, use_site = use_site_now)
    if (!is.null(m2)) {
      res_rob[[i]] <- m2 %>%
        mutate(
          table = table_name,
          hemi = meta$hemi[1],
          atlas = meta$atlas[1],
          roi = r,
          n_obs = nrow(dat_roi),
          n_sub = n_distinct(dat_roi$ident),
          n_ge2 = {
            tmp <- dat_roi %>% distinct(ident, wave) %>% count(ident, name = "n_waves")
            sum(tmp$n_waves >= 2)
          },
          n_ge3 = {
            tmp <- dat_roi %>% distinct(ident, wave) %>% count(ident, name = "n_waves")
            sum(tmp$n_waves >= 3)
          }
        )
    }
  }
  
  main_tbl_raw <- bind_rows(res_main)
  rob_tbl_raw  <- bind_rows(res_rob)
  
  write_csv(main_tbl_raw, file.path(OUT_TMP, paste0("MAIN_RAW_", table_name, ".csv")))
  write_csv(rob_tbl_raw,  file.path(OUT_TMP, paste0("ROB_RAW_", table_name, ".csv")))
  write_csv(sample_qc,    file.path(OUT_TMP, paste0("SAMPLE_QC_", table_name, ".csv")))
  
  list(
    table_name = table_name,
    wide_m = wide_m,
    sample_qc = sample_qc,
    main_raw = main_tbl_raw,
    rob_raw = rob_tbl_raw
  )
}

# ==============================================================================
# 5) RUN LH THEN RH
# ==============================================================================
msg("\n[4/10] Running LH ...")
res_lh <- run_one_hemi(file.path(PATH_FS_TABLES, "lh_aparc_thickness.tsv"))

msg("\n[5/10] Running RH ...")
res_rh <- run_one_hemi(file.path(PATH_FS_TABLES, "rh_aparc_thickness.tsv"))

# ==============================================================================
# 6) MERGE LH + RH
# ==============================================================================
msg("\n[6/10] Merging LH + RH ...")

MAIN_RES <- bind_rows(res_lh$main_raw, res_rh$main_raw) %>%
  filter(!is.na(term), !is.na(estimate), !is.na(p.value))

ROB_RES <- bind_rows(res_lh$rob_raw, res_rh$rob_raw) %>%
  filter(!is.na(term), !is.na(estimate), !is.na(p.value))

SAMPLE_QC_ALL <- bind_rows(res_lh$sample_qc, res_rh$sample_qc)

if (nrow(MAIN_RES) == 0) stop("MAIN_RES is empty after merging LH and RH.")
if (nrow(ROB_RES) == 0) stop("ROB_RES is empty after merging LH and RH.")

write_csv(MAIN_RES, file.path(OUT_MAIN, "01_MAIN_all_terms_raw_LH_RH.csv"))
write_csv(ROB_RES,  file.path(OUT_ROB,  "01_ROBUST_all_terms_raw_LH_RH.csv"))
write_csv(SAMPLE_QC_ALL, file.path(OUT_QC, "00_QC_ROI_SAMPLE_OVERVIEW.csv"))

MAIN_RES <- MAIN_RES %>%
  group_by(hemi, term) %>%
  mutate(pFDR_hemi = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(pFDR_global = p.adjust(p.value, method = "fdr")) %>%
  ungroup()

ROB_RES <- ROB_RES %>%
  group_by(hemi, term) %>%
  mutate(pFDR_hemi = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  group_by(term) %>%
  mutate(pFDR_global = p.adjust(p.value, method = "fdr")) %>%
  ungroup()

write_csv(MAIN_RES, file.path(OUT_MAIN, "02_MAIN_all_terms_FDR_LH_RH.csv"))
write_csv(ROB_RES,  file.path(OUT_ROB,  "02_ROBUST_all_terms_FDR_LH_RH.csv"))

# ==============================================================================
# 7) TABLES / SHORTLISTS
# ==============================================================================
msg("\n[7/10] Building manuscript tables and shortlists ...")

MAIN_KEY <- MAIN_RES %>%
  filter(term %in% KEY_TERMS_MAIN) %>%
  mutate(
    term_pretty = term_label_main(term),
    beta = estimate,
    se = std.error,
    t = statistic,
    p = p.value,
    p_fdr_hemi = pFDR_hemi,
    p_fdr_global = pFDR_global,
    roi_label = pretty_roi_label(roi)
  ) %>%
  select(
    table, hemi, atlas, roi, roi_label, term, term_pretty,
    beta, se, conf.low, conf.high, t, p, p_fdr_hemi, p_fdr_global,
    model, model_random, singular, converged, n_sub, n_obs, n_ge2, n_ge3
  ) %>%
  arrange(term, p_fdr_global, desc(abs(beta))) %>%
  distinct(table, roi, term, .keep_all = TRUE)

ROB_KEY <- ROB_RES %>%
  filter(term %in% KEY_TERMS_ROB) %>%
  mutate(
    term_pretty = term_label_rob(term),
    beta = estimate,
    se = std.error,
    t = statistic,
    p = p.value,
    p_fdr_hemi = pFDR_hemi,
    p_fdr_global = pFDR_global,
    roi_label = pretty_roi_label(roi)
  ) %>%
  select(
    table, hemi, atlas, roi, roi_label, term, term_pretty,
    beta, se, conf.low, conf.high, t, p, p_fdr_hemi, p_fdr_global,
    model, model_random, singular, converged, n_sub, n_obs, n_ge2, n_ge3
  ) %>%
  arrange(term, p_fdr_global, desc(abs(beta))) %>%
  distinct(table, roi, term, .keep_all = TRUE)

write_csv(MAIN_KEY, file.path(OUT_TAB, "01_MAIN_MANUSCRIPT_KEY.csv"))
write_csv(ROB_KEY,  file.path(OUT_TAB, "02_ROBUST_MANUSCRIPT_KEY.csv"))

make_shortlist <- function(df, term, fdr_col = "p_fdr_global",
                           fdr_main = 0.05, fdr_fallback = 0.10){
  out <- df %>%
    filter(
      term == !!term,
      is_cortical_roi(roi),
      !is.na(beta),
      !is.na(.data[[fdr_col]]),
      .data[[fdr_col]] < fdr_main
    )
  if (nrow(out) == 0) {
    out <- df %>%
      filter(
        term == !!term,
        is_cortical_roi(roi),
        !is.na(beta),
        !is.na(.data[[fdr_col]]),
        .data[[fdr_col]] < fdr_fallback
      )
  }
  out %>% arrange(.data[[fdr_col]], desc(abs(beta)))
}

SHORT_cPES_GLOBAL     <- make_shortlist(MAIN_KEY, TERM_MAIN_cPES,     "p_fdr_global", FDR_MAIN, FDR_FALLBACK)
SHORT_AXE_GLOBAL     <- make_shortlist(MAIN_KEY, TERM_MAIN_AXE,     "p_fdr_global", FDR_MAIN, FDR_FALLBACK)
SHORT_AXE_SEX_GLOBAL <- make_shortlist(MAIN_KEY, TERM_MAIN_AXE_SEX, "p_fdr_global", FDR_MAIN, FDR_FALLBACK)

SHORT_cPES_HEMI     <- make_shortlist(MAIN_KEY, TERM_MAIN_cPES,     "p_fdr_hemi", FDR_MAIN, FDR_FALLBACK)
SHORT_AXE_HEMI     <- make_shortlist(MAIN_KEY, TERM_MAIN_AXE,     "p_fdr_hemi", FDR_MAIN, FDR_FALLBACK)
SHORT_AXE_SEX_HEMI <- make_shortlist(MAIN_KEY, TERM_MAIN_AXE_SEX, "p_fdr_hemi", FDR_MAIN, FDR_FALLBACK)

write_csv(SHORT_cPES_GLOBAL,     file.path(OUT_TAB, "03_SHORTLIST_MAIN_cPES_GLOBAL.csv"))
write_csv(SHORT_AXE_GLOBAL,     file.path(OUT_TAB, "04_SHORTLIST_AGE_X_cPES_GLOBAL.csv"))
write_csv(SHORT_AXE_SEX_GLOBAL, file.path(OUT_TAB, "05_SHORTLIST_AGE_X_cPES_X_SEX_GLOBAL.csv"))

write_csv(SHORT_cPES_HEMI,     file.path(OUT_TAB, "06_SHORTLIST_MAIN_cPES_HEMI.csv"))
write_csv(SHORT_AXE_HEMI,     file.path(OUT_TAB, "07_SHORTLIST_AGE_X_cPES_HEMI.csv"))
write_csv(SHORT_AXE_SEX_HEMI, file.path(OUT_TAB, "08_SHORTLIST_AGE_X_cPES_X_SEX_HEMI.csv"))

ROBUST_COMPARE_AXE <- SHORT_AXE_HEMI %>%
  select(table, roi, roi_label, beta_main = beta, p_fdr_main_hemi = p_fdr_hemi) %>%
  left_join(
    ROB_KEY %>%
      filter(term == TERM_ROB_AXE) %>%
      select(table, roi, beta_rob = beta, p_fdr_rob_global = p_fdr_global),
    by = c("table", "roi")
  ) %>%
  arrange(p_fdr_main_hemi)

write_csv(ROBUST_COMPARE_AXE, file.path(OUT_TAB, "09_MAIN_vs_ROBUST_COMPARE_AGE_X_cPES_HEMI_PRIMARY.csv"))

# PRIMARY FIGURES NOW USE HEMISPHERE-LEVEL FDR
SHORT_AXE <- SHORT_AXE_HEMI
SHORT_AXE_SEX <- SHORT_AXE_SEX_HEMI

# ==============================================================================
# 8) BRAIN MAP
# ==============================================================================
msg("\n[8/10] Building brain map ...")

BRAINMAP_AXE <- MAIN_KEY %>%
  filter(term == TERM_MAIN_AXE, is_cortical_roi(roi)) %>%
  filter(!str_detect(tolower(roi), "meanthickness")) %>%
  mutate(
    hemi_fs = ifelse(str_starts(roi, "lh_"), "lh", "rh"),
    hemi_ggseg = get_hemi_from_roi(roi),
    region_fs = get_region_from_roi(roi),
    region_ggseg = to_ggseg_region(roi),
    significant_fdr05_global = p_fdr_global < 0.05,
    significant_fdr10_global = p_fdr_global < 0.10,
    significant_fdr05_hemi = p_fdr_hemi < 0.05,
    significant_fdr10_hemi = p_fdr_hemi < 0.10,
    t_plot_global = ifelse(significant_fdr05_global, t, NA_real_),
    t_plot_hemi = ifelse(significant_fdr05_hemi, t, NA_real_),
    beta_plot = ifelse(significant_fdr05_hemi, beta, NA_real_)
  ) %>%
  arrange(p_fdr_hemi)

write_csv(BRAINMAP_AXE, file.path(OUT_MAPS, "01_brainmap_age_x_cpes_values.csv"))

if (HAS_GGSEG && nrow(BRAINMAP_AXE) > 0) {
  msg("  ggseg detected -> generating primary hemisphere-level brain map ...")
  
  suppressMessages({
    library(ggseg)
  })
  
  map_df <- BRAINMAP_AXE %>%
    transmute(
      region = region_ggseg,
      hemi = hemi_ggseg,
      fill_value = t_plot_hemi,
      beta = beta,
      p_fdr_global = p_fdr_global,
      p_fdr_hemi = p_fdr_hemi
    ) %>%
    filter(region != "banksts") %>%
    distinct()
  
  max_abs_t <- max(abs(BRAINMAP_AXE$t[is.finite(BRAINMAP_AXE$t)]), na.rm = TRUE)
  if (!is.finite(max_abs_t) || is.na(max_abs_t)) max_abs_t <- 3
  
  p_map <- ggplot(map_df) +
    ggseg::geom_brain(
      atlas = ggseg::dk(),
      aes(fill = fill_value),
      color = "grey35",
      size = 0.15
    ) +
    scale_fill_gradient2(
      low = COL_LOW,
      mid = "grey95",
      high = COL_HIGH,
      midpoint = 0,
      na.value = "grey85",
      limits = c(-max_abs_t, max_abs_t),
      name = "t-value\n(age × cPES)"
    ) +
    labs(
      title = "Cortical regions showing age × cPES effects",
      subtitle = "Primary map uses hemisphere-level FDR-significant ROIs; non-significant regions are shown in grey"
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave(
    file.path(OUT_FIGS, "FIG_A_brainmap_age_x_cpes_unified.png"),
    p_map, width = 10, height = 8, dpi = DPI_FIG, bg = "white"
  )
} else {
  msg("  ggseg not available or no age × cPES hits -> exported CSV only.")
}

# ==============================================================================
# 9) TRAJECTORIES + DELTA SLOPE FOR ALL SHORTLISTED ROIS
# ==============================================================================
msg("\n[9/10] Building trajectories and delta slope for all shortlisted ROIs ...")

compute_conditional_age_slope <- function(model_obj, cpes_values, sex = "Male"){
  tb <- tidy(model_obj, effects = "fixed", conf.int = FALSE)
  
  get_coef <- function(term){
    val <- tb$estimate[tb$term == term][1]
    if (length(val) == 0 || is.na(val)) 0 else val
  }
  
  vc <- as.matrix(vcov(model_obj))
  rn <- rownames(vc)
  
  b_age      <- get_coef("age_wave")
  b_axe      <- get_coef("age_wave:cPES_z")
  b_age_fem  <- get_coef("age_wave:genderFemale")
  b_axesex   <- get_coef("age_wave:cPES_z:genderFemale")
  
  female_ind <- ifelse(sex == "Female", 1, 0)
  
  out <- lapply(cpes_values, function(e){
    
    slope <- b_age + b_axe * e + b_age_fem * female_ind + b_axesex * e * female_ind
    
    coef_vec <- rep(0, length(rn))
    names(coef_vec) <- rn
    if ("age_wave" %in% rn) coef_vec["age_wave"] <- 1
    if ("age_wave:cPES_z" %in% rn) coef_vec["age_wave:cPES_z"] <- e
    if (female_ind == 1 && "age_wave:genderFemale" %in% rn) coef_vec["age_wave:genderFemale"] <- 1
    if (female_ind == 1 && "age_wave:cPES_z:genderFemale" %in% rn) coef_vec["age_wave:cPES_z:genderFemale"] <- e
    
    se <- sqrt(as.numeric(t(coef_vec) %*% vc %*% coef_vec))
    
    tibble(
      cPES_z = e,
      sex = sex,
      slope_age = slope,
      conf.low = slope - 1.96 * se,
      conf.high = slope + 1.96 * se
    )
  })
  
  bind_rows(out)
}

if (nrow(SHORT_AXE) > 0) {
  
  write_csv(SHORT_AXE, file.path(OUT_TAB, "10_ALL_SHORTLISTED_AGE_X_cPES_ROIS_FOR_FIGURES.csv"))
  
  for (i in seq_len(nrow(SHORT_AXE))) {
    
    TABLE_PICK <- SHORT_AXE$table[i]
    ROI_PICK   <- SHORT_AXE$roi[i]
    ROI_LABEL  <- SHORT_AXE$roi_label[i]
    
    wide_pick <- read_csv(
      file.path(OUT_QC, paste0("00_", TABLE_PICK, "_wide_merged.csv")),
      show_col_types = FALSE
    ) %>% clean_names()
    
    roi_col <- resolve_roi_col(wide_pick, ROI_PICK)
    cpes_col <- get_cpes_col_after_clean(wide_pick)
    if (is.na(cpes_col)) stop("cPES_z / cpes_z column not found after clean_names() in: ", TABLE_PICK)
    
    dat_pick <- wide_pick %>%
      transmute(
        ident = as.character(ident),
        age_wave = as.numeric(age_wave),
        cPES_z = as.numeric(.data[[cpes_col]]),
        gender = factor(gender, levels = c("Male", "Female")),
        site_original = factor(site_original),
        value = suppressWarnings(as.numeric(.data[[roi_col]]))
      ) %>%
      filter(!is.na(value), !is.na(age_wave), !is.na(cPES_z), !is.na(gender))
    
    use_site_pick <- USE_SITE && (nlevels(dat_pick$site_original) >= 2)
    ref_site <- if (use_site_pick) pick_reference_site(dat_pick$site_original) else NA_character_
    
    f_plot <- as.formula(paste0(
      "value ~ age_wave * cPES_z * gender",
      if (use_site_pick) " + site_original" else "",
      " + (1 | ident)"
    ))
    
    fit_plot <- lmer(f_plot, data = dat_pick, REML = FALSE)
    
    amin <- min(dat_pick$age_wave, na.rm = TRUE)
    amax <- max(dat_pick$age_wave, na.rm = TRUE)
    
    nd <- expand.grid(
      age_wave = seq(amin, amax, length.out = 120),
      cPES_z = c(-1, 0, 1),
      gender = factor("Male", levels = c("Male", "Female")),
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    
    if (use_site_pick) {
      nd$site_original <- factor(ref_site, levels = levels(dat_pick$site_original))
    }
    
    pred_link <- predict(fit_plot, newdata = nd, re.form = NA, se.fit = TRUE, allow.new.levels = TRUE)
    nd$predicted <- as.numeric(pred_link$fit)
    nd$se_fit <- as.numeric(pred_link$se.fit)
    nd$conf.low <- nd$predicted - 1.96 * nd$se_fit
    nd$conf.high <- nd$predicted + 1.96 * nd$se_fit
    
    nd <- nd %>%
      mutate(
        cPES_group = factor(
          cPES_z,
          levels = c(-1, 0, 1),
          labels = c("Low cPES (-1 SD)", "Mean cPES", "High cPES (+1 SD)")
        ),
        roi_label = ROI_LABEL
      )
    
    plot_lowhigh <- nd %>% filter(cPES_group %in% c("Low cPES (-1 SD)", "High cPES (+1 SD)"))
    
    n_low <- dat_pick %>%
      filter(cPES_z <= -1) %>%
      distinct(ident) %>%
      nrow()
    
    n_high <- dat_pick %>%
      filter(cPES_z >= 1) %>%
      distinct(ident) %>%
      nrow()
    
    subtitle_n <- paste0("Low cPES (n=", n_low, ") | High cPES (n=", n_high, ")")
    
    legend_labs <- c(
      paste0("Low cPES (-1 SD), n=", n_low),
      paste0("High cPES (+1 SD), n=", n_high)
    )
    
    plot_lowhigh <- plot_lowhigh %>%
      mutate(
        cPES_group_n = factor(
          cPES_group,
          levels = c("Low cPES (-1 SD)", "High cPES (+1 SD)"),
          labels = legend_labs
        )
      )
    
    dat_lowhigh_points <- dat_pick %>%
      mutate(
        cPES_group = case_when(
          cPES_z <= -1 ~ "Low cPES (-1 SD)",
          cPES_z >=  1 ~ "High cPES (+1 SD)",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(cPES_group)) %>%
      mutate(
        cPES_group = factor(
          cPES_group,
          levels = c("Low cPES (-1 SD)", "High cPES (+1 SD)"),
          labels = legend_labs
        ),
        age_wave_j = jitter(age_wave, amount = 0.06),
        value_j = jitter(value, amount = 0)
      )
    
    # FIG_B1 — low/high only, N in subtitle
    p_traj_b1 <- ggplot() +
      geom_ribbon(
        data = plot_lowhigh,
        aes(x = age_wave, ymin = conf.low, ymax = conf.high, fill = cPES_group),
        alpha = 0.18
      ) +
      geom_line(
        data = plot_lowhigh,
        aes(x = age_wave, y = predicted, color = cPES_group),
        linewidth = 1.25
      ) +
      theme_paper(base_size = 13) +
      labs(
        x = "Age (years)",
        y = "Cortical thickness",
        title = ROI_LABEL,
        subtitle = subtitle_n
      ) +
      scale_color_manual(
        name = NULL,
        values = c(
          "Low cPES (-1 SD)" = COL_LOW,
          "High cPES (+1 SD)" = COL_HIGH
        )
      ) +
      scale_fill_manual(
        name = NULL,
        values = c(
          "Low cPES (-1 SD)" = COL_LOW,
          "High cPES (+1 SD)" = COL_HIGH
        )
      ) +
      theme(
        legend.position = "bottom",
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(
      file.path(OUT_FIGS, paste0("FIG_B1_TRAJ_LOWHIGH_SUBTITLE_N_", i, "_", TABLE_PICK, "_", ROI_PICK, ".png")),
      p_traj_b1, width = 7.0, height = 5.4, dpi = DPI_FIG, bg = "white"
    )
    
    # FIG_B2 — low/high only, N in legend
    p_traj_b2 <- ggplot() +
      geom_ribbon(
        data = plot_lowhigh,
        aes(x = age_wave, ymin = conf.low, ymax = conf.high, fill = cPES_group_n),
        alpha = 0.18
      ) +
      geom_line(
        data = plot_lowhigh,
        aes(x = age_wave, y = predicted, color = cPES_group_n),
        linewidth = 1.25
      ) +
      theme_paper(base_size = 13) +
      labs(
        x = "Age (years)",
        y = "Cortical thickness",
        title = ROI_LABEL,
        subtitle = "Predicted trajectories from the primary LMM"
      ) +
      scale_color_manual(
        name = NULL,
        values = setNames(c(COL_LOW, COL_HIGH), legend_labs)
      ) +
      scale_fill_manual(
        name = NULL,
        values = setNames(c(COL_LOW, COL_HIGH), legend_labs)
      ) +
      theme(
        legend.position = "bottom",
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(
      file.path(OUT_FIGS, paste0("FIG_B2_TRAJ_LOWHIGH_LEGEND_N_", i, "_", TABLE_PICK, "_", ROI_PICK, ".png")),
      p_traj_b2, width = 7.0, height = 5.4, dpi = DPI_FIG, bg = "white"
    )
    
    # FIG_B3 — low/high only, observed points + model trajectories + N in legend
    p_traj_b3 <- ggplot() +
      geom_point(
        data = dat_lowhigh_points,
        aes(x = age_wave_j, y = value_j, color = cPES_group),
        alpha = 0.28,
        size = 1.2
      ) +
      geom_ribbon(
        data = plot_lowhigh,
        aes(x = age_wave, ymin = conf.low, ymax = conf.high, fill = cPES_group_n),
        alpha = 0.16
      ) +
      geom_line(
        data = plot_lowhigh,
        aes(x = age_wave, y = predicted, color = cPES_group_n),
        linewidth = 1.25
      ) +
      theme_paper(base_size = 13) +
      labs(
        x = "Age (years)",
        y = "Cortical thickness",
        title = ROI_LABEL,
        subtitle = "Observed phenotype groups overlaid with model-predicted trajectories"
      ) +
      scale_color_manual(
        name = NULL,
        values = setNames(c(COL_LOW, COL_HIGH), legend_labs),
        breaks = legend_labs
      ) +
      scale_fill_manual(
        name = NULL,
        values = setNames(c(COL_LOW, COL_HIGH), legend_labs),
        breaks = legend_labs
      ) +
      theme(
        legend.position = "bottom",
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(
      file.path(OUT_FIGS, paste0("FIG_B3_TRAJ_LOWHIGH_POINTS_", i, "_", TABLE_PICK, "_", ROI_PICK, ".png")),
      p_traj_b3, width = 7.0, height = 5.4, dpi = DPI_FIG, bg = "white"
    )
    
    cpes_seq <- seq(
      quantile(dat_pick$cPES_z, 0.05, na.rm = TRUE),
      quantile(dat_pick$cPES_z, 0.95, na.rm = TRUE),
      length.out = 120
    )
    
    slope_male <- compute_conditional_age_slope(fit_plot, cpes_seq, sex = "Male")
    slope_female <- compute_conditional_age_slope(fit_plot, cpes_seq, sex = "Female")
    slope_df <- bind_rows(slope_male, slope_female)
    
    write_csv(
      slope_df,
      file.path(OUT_TAB, paste0("11_DELTA_SLOPE_DATA_ALL_", i, "_", TABLE_PICK, "_", ROI_PICK, ".csv"))
    )
    
    p_delta <- slope_male %>%
      ggplot(aes(x = cPES_z, y = slope_age)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.18, fill = "grey70") +
      geom_line(linewidth = 1.1, color = "black") +
      geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
      theme_paper(base_size = 12) +
      labs(
        x = "Exposome Risk Score (cPES z)",
        y = "Conditional age slope on thickness",
        title = paste0("Conditional age effect across cPES: ", ROI_LABEL),
        subtitle = "More negative slopes indicate steeper age-related cortical thinning"
      )
    
    ggsave(
      file.path(OUT_FIGS, paste0("FIG_C_DELTA_SLOPE_ALL_", i, "_", TABLE_PICK, "_", ROI_PICK, ".png")),
      p_delta, width = 7.0, height = 5.3, dpi = DPI_FIG, bg = "white"
    )
  }
  
} else {
  msg("  No age × cPES cortical ROIs -> skipping trajectories/delta slope.")
}

# Sex-specific delta slopes for all shortlisted age x cPES x sex ROIs
if (nrow(SHORT_AXE_SEX) > 0) {
  
  write_csv(SHORT_AXE_SEX, file.path(OUT_TAB, "12_ALL_SHORTLISTED_AGE_X_cPES_X_SEX_ROIS_FOR_FIGURES.csv"))
  
  for (i in seq_len(nrow(SHORT_AXE_SEX))) {
    
    TABLE_PICK <- SHORT_AXE_SEX$table[i]
    ROI_PICK   <- SHORT_AXE_SEX$roi[i]
    ROI_LABEL  <- SHORT_AXE_SEX$roi_label[i]
    
    wide_pick <- read_csv(
      file.path(OUT_QC, paste0("00_", TABLE_PICK, "_wide_merged.csv")),
      show_col_types = FALSE
    ) %>% clean_names()
    
    roi_col <- resolve_roi_col(wide_pick, ROI_PICK)
    cpes_col <- get_cpes_col_after_clean(wide_pick)
    if (is.na(cpes_col)) stop("cPES_z / cpes_z column not found after clean_names() in: ", TABLE_PICK)
    
    dat_pick <- wide_pick %>%
      transmute(
        ident = as.character(ident),
        age_wave = as.numeric(age_wave),
        cPES_z = as.numeric(.data[[cpes_col]]),
        gender = factor(gender, levels = c("Male", "Female")),
        site_original = factor(site_original),
        value = suppressWarnings(as.numeric(.data[[roi_col]]))
      ) %>%
      filter(!is.na(value), !is.na(age_wave), !is.na(cPES_z), !is.na(gender))
    
    use_site_pick <- USE_SITE && (nlevels(dat_pick$site_original) >= 2)
    
    f_plot <- as.formula(paste0(
      "value ~ age_wave * cPES_z * gender",
      if (use_site_pick) " + site_original" else "",
      " + (1 | ident)"
    ))
    
    fit_plot <- lmer(f_plot, data = dat_pick, REML = FALSE)
    
    cpes_seq <- seq(
      quantile(dat_pick$cPES_z, 0.05, na.rm = TRUE),
      quantile(dat_pick$cPES_z, 0.95, na.rm = TRUE),
      length.out = 120
    )
    
    slope_male <- compute_conditional_age_slope(fit_plot, cpes_seq, sex = "Male")
    slope_female <- compute_conditional_age_slope(fit_plot, cpes_seq, sex = "Female")
    slope_df <- bind_rows(slope_male, slope_female)
    
    p_delta_sex <- slope_df %>%
      ggplot(aes(x = cPES_z, y = slope_age, linetype = sex, fill = sex)) +
      geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.12, colour = NA) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
      theme_paper(base_size = 12) +
      labs(
        x = "Exposome Risk Score (cPES z)",
        y = "Conditional age slope on thickness",
        title = paste0("Sex-specific conditional age slopes: ", ROI_LABEL),
        subtitle = "Supplementary figure for ROIs with significant age × cPES × sex interaction",
        linetype = "",
        fill = ""
      ) +
      theme(legend.position = "bottom")
    
    ggsave(
      file.path(OUT_FIGS, paste0("SUPP_DELTA_SLOPE_SEX_ALL_", i, "_", TABLE_PICK, "_", ROI_PICK, ".png")),
      p_delta_sex, width = 7.3, height = 5.6, dpi = DPI_FIG, bg = "white"
    )
  }
}

# ==============================================================================
# 10) SUPPLEMENTARY FIGURES + DEVELOPMENTAL YEARS
# ==============================================================================
msg("\n[10/10] Supplementary outputs ...")

if (nrow(SHORT_AXE) > 0) {
  
  forest_df <- MAIN_KEY %>%
    filter(term == TERM_MAIN_AXE, is_cortical_roi(roi), !is.na(p_fdr_hemi), p_fdr_hemi < max(FDR_MAIN, FDR_FALLBACK)) %>%
    arrange(p_fdr_hemi, desc(abs(beta))) %>%
    slice_head(n = TOP_PER_TERM_FOREST) %>%
    mutate(
      label = paste0(ifelse(str_starts(roi, "lh_"), "LH", "RH"), " | ", roi_label),
      ci_low = beta - 1.96 * se,
      ci_high = beta + 1.96 * se
    )
  
  p_forest <- ggplot(forest_df, aes(x = beta, y = fct_reorder(label, beta))) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_errorbar(aes(y = fct_reorder(label, beta), xmin = ci_low, xmax = ci_high),
                  orientation = "y", width = 0.2) +
    geom_point(size = 2) +
    theme_paper(base_size = 11) +
    labs(
      x = "Beta (approx. 95% CI)",
      y = "",
      title = "Supplementary: top cortical regions showing age × cPES effects"
    )
  
  ggsave(
    file.path(OUT_FIGS, "SUPP_FIG_1_FOREST_AGE_X_cPES.png"),
    p_forest, width = 9.5, height = 7.5, dpi = DPI_FIG, bg = "white"
  )
}

compute_age_equivalent_one <- function(dat_roi, use_site = TRUE){
  
  f <- as.formula(
    paste0(
      "value ~ age_wave * cPES_z * gender",
      if (use_site) " + site_original" else "",
      " + (1 | ident)"
    )
  )
  
  m <- tryCatch(lmer(f, data = dat_roi, REML = FALSE), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  
  tb <- tidy(m, effects = "fixed", conf.int = TRUE)
  
  beta_age <- tb$estimate[tb$term == "age_wave"][1]
  beta_axe <- tb$estimate[tb$term == "age_wave:cPES_z"][1]
  
  if (length(beta_age) == 0 || length(beta_axe) == 0 ||
      is.na(beta_age) || is.na(beta_axe)) return(NULL)
  
  q25 <- as.numeric(quantile(dat_roi$cPES_z, 0.25, na.rm = TRUE))
  q75 <- as.numeric(quantile(dat_roi$cPES_z, 0.75, na.rm = TRUE))
  delta_cpes <- q75 - q25
  
  extra_thinning_window <- beta_axe * delta_cpes * 5
  age_equiv_years <- ifelse(
    abs(beta_age) > 1e-10,
    extra_thinning_window / abs(beta_age),
    NA_real_
  )
  
  tibble(
    beta_age = beta_age,
    beta_ageXcPES = beta_axe,
    cpes_q25 = q25,
    cpes_q75 = q75,
    delta_cpes = delta_cpes,
    age_window_years = 5,
    extra_thinning_over_window = extra_thinning_window,
    developmental_years_equivalent = age_equiv_years,
    interpretation_note = "Heuristic supplementary translation only; not for primary claim."
  )
}

impact_list <- list()

if (nrow(SHORT_AXE) > 0) {
  for (i in seq_len(nrow(SHORT_AXE))) {
    tb <- SHORT_AXE$table[i]
    rr <- SHORT_AXE$roi[i]
    
    wide_tmp <- read_csv(
      file.path(OUT_QC, paste0("00_", tb, "_wide_merged.csv")),
      show_col_types = FALSE
    ) %>% clean_names()
    
    roi_col <- resolve_roi_col(wide_tmp, rr)
    cpes_col <- get_cpes_col_after_clean(wide_tmp)
    if (is.na(cpes_col)) stop("cPES_z / cpes_z column not found after clean_names() in: ", tb)
    
    dat_roi <- wide_tmp %>%
      transmute(
        ident = as.character(ident),
        age_wave = as.numeric(age_wave),
        cPES_z = as.numeric(.data[[cpes_col]]),
        gender = factor(gender, levels = c("Male", "Female")),
        site_original = factor(site_original),
        value = suppressWarnings(as.numeric(.data[[roi_col]]))
      ) %>%
      filter(!is.na(value), !is.na(age_wave), !is.na(cPES_z), !is.na(gender))
    
    out <- compute_age_equivalent_one(dat_roi, use_site = USE_SITE)
    
    if (!is.null(out)) {
      impact_list[[paste(tb, rr, sep = "__")]] <- out %>%
        mutate(table = tb, roi = rr, roi_label = pretty_roi_label(rr))
    }
  }
}

IMPACT_EQ <- bind_rows(impact_list) %>%
  left_join(
    SHORT_AXE %>% select(table, roi, term, term_pretty, beta, se, p_fdr_hemi, p_fdr_global),
    by = c("table", "roi")
  ) %>%
  arrange(p_fdr_hemi)

write_csv(IMPACT_EQ, file.path(OUT_TAB, "13_developmental_years_equivalent_SUPPLEMENTARY.csv"))

msg("\nDONE. Outputs saved in: %s", normalizePath(OUT_BASE))

################################################################################
# EXTRA BLOCK — HEMISPHERIC BRAIN MAP + AGE-SPECIFIC DIFFERENCE PLOTS
################################################################################

# ============================================================
# A) Brain map using hemisphere-level FDR
# ============================================================
BRAINMAP_AXE_HEMI <- MAIN_KEY %>%
  filter(term == TERM_MAIN_AXE, is_cortical_roi(roi)) %>%
  filter(!str_detect(tolower(roi), "meanthickness")) %>%
  mutate(
    hemi_fs = ifelse(str_starts(roi, "lh_"), "lh", "rh"),
    hemi_ggseg = get_hemi_from_roi(roi),
    region_fs = get_region_from_roi(roi),
    region_ggseg = to_ggseg_region(roi),
    significant_fdr05_hemi = p_fdr_hemi < 0.05,
    t_plot = ifelse(significant_fdr05_hemi, t, NA_real_)
  ) %>%
  filter(region_ggseg != "banksts") %>%
  arrange(p_fdr_hemi)

write_csv(BRAINMAP_AXE_HEMI, file.path(OUT_MAPS, "02_brainmap_age_x_cpes_values_HEMI_FDR.csv"))

if (HAS_GGSEG && nrow(BRAINMAP_AXE_HEMI) > 0) {
  suppressMessages(library(ggseg))
  
  map_df_hemi <- BRAINMAP_AXE_HEMI %>%
    transmute(
      region = region_ggseg,
      hemi = hemi_ggseg,
      fill_value = t_plot
    ) %>%
    distinct()
  
  max_abs_t <- max(abs(BRAINMAP_AXE_HEMI$t[is.finite(BRAINMAP_AXE_HEMI$t)]), na.rm = TRUE)
  if (!is.finite(max_abs_t) || is.na(max_abs_t)) max_abs_t <- 3
  
  p_map_hemi <- ggplot(map_df_hemi) +
    ggseg::geom_brain(
      atlas = ggseg::dk(),
      aes(fill = fill_value),
      color = "grey35",
      size = 0.15
    ) +
    scale_fill_gradient2(
      low = COL_LOW,
      mid = "grey95",
      high = COL_HIGH,
      midpoint = 0,
      na.value = "grey85",
      limits = c(-max_abs_t, max_abs_t),
      name = "t-value\n(age × cPES)"
    ) +
    labs(
      title = "Cortical regions showing age × cPES effects",
      subtitle = "Map shows hemisphere-level FDR-significant ROIs; non-significant regions are shown in grey"
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  ggsave(
    file.path(OUT_FIGS, "FIG_A_brainmap_age_x_cpes_HEMI_FDR.png"),
    p_map_hemi, width = 10, height = 8, dpi = DPI_FIG, bg = "white"
  )
}

# ============================================================
# B) Age-specific difference plot (High cPES - Low cPES)
# ============================================================
for (i in seq_len(nrow(SHORT_AXE))) {
  
  TABLE_PICK <- SHORT_AXE$table[i]
  ROI_PICK   <- SHORT_AXE$roi[i]
  ROI_LABEL  <- SHORT_AXE$roi_label[i]
  
  wide_pick <- read_csv(
    file.path(OUT_QC, paste0("00_", TABLE_PICK, "_wide_merged.csv")),
    show_col_types = FALSE
  ) %>% clean_names()
  
  roi_col <- resolve_roi_col(wide_pick, ROI_PICK)
  cpes_col <- get_cpes_col_after_clean(wide_pick)
  
  dat_pick <- wide_pick %>%
    transmute(
      ident = as.character(ident),
      age_wave = as.numeric(age_wave),
      cPES_z = as.numeric(.data[[cpes_col]]),
      gender = factor(gender, levels = c("Male", "Female")),
      site_original = factor(site_original),
      value = suppressWarnings(as.numeric(.data[[roi_col]]))
    ) %>%
    filter(!is.na(value), !is.na(age_wave), !is.na(cPES_z), !is.na(gender))
  
  use_site_pick <- USE_SITE && (nlevels(dat_pick$site_original) >= 2)
  ref_site <- if (use_site_pick) pick_reference_site(dat_pick$site_original) else NA_character_
  
  f_plot <- as.formula(paste0(
    "value ~ age_wave * cPES_z * gender",
    if (use_site_pick) " + site_original" else "",
    " + (1 | ident)"
  ))
  
  fit_plot <- lmer(f_plot, data = dat_pick, REML = FALSE)
  
  amin <- min(dat_pick$age_wave, na.rm = TRUE)
  amax <- max(dat_pick$age_wave, na.rm = TRUE)
  
  nd_low <- data.frame(
    age_wave = seq(amin, amax, length.out = 120),
    cPES_z = -1,
    gender = factor("Male", levels = c("Male", "Female"))
  )
  
  nd_high <- data.frame(
    age_wave = seq(amin, amax, length.out = 120),
    cPES_z = 1,
    gender = factor("Male", levels = c("Male", "Female"))
  )
  
  if (use_site_pick) {
    nd_low$site_original  <- factor(ref_site, levels = levels(dat_pick$site_original))
    nd_high$site_original <- factor(ref_site, levels = levels(dat_pick$site_original))
  }
  
  pred_low  <- predict(fit_plot, newdata = nd_low,  re.form = NA, se.fit = TRUE, allow.new.levels = TRUE)
  pred_high <- predict(fit_plot, newdata = nd_high, re.form = NA, se.fit = TRUE, allow.new.levels = TRUE)
  
  diff_df <- tibble(
    age_wave = nd_low$age_wave,
    diff_high_minus_low = as.numeric(pred_high$fit) - as.numeric(pred_low$fit),
    se_diff = sqrt(as.numeric(pred_high$se.fit)^2 + as.numeric(pred_low$se.fit)^2),
    conf.low = diff_high_minus_low - 1.96 * se_diff,
    conf.high = diff_high_minus_low + 1.96 * se_diff,
    roi_label = ROI_LABEL
  )
  
  write_csv(
    diff_df,
    file.path(OUT_TAB, paste0("14_AGE_SPECIFIC_DIFF_", i, "_", TABLE_PICK, "_", ROI_PICK, ".csv"))
  )
  
  p_diff <- ggplot(diff_df, aes(x = age_wave, y = diff_high_minus_low)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.18, fill = "grey75") +
    geom_line(linewidth = 1.15, color = COL_HIGH) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    theme_paper(base_size = 12) +
    labs(
      x = "Age (years)",
      y = "Predicted thickness difference\n(High cPES - Low cPES)",
      title = paste0("Age-specific cPES difference: ", ROI_LABEL),
      subtitle = "Negative values indicate lower thickness in the high-cPES trajectory"
    )
  
  ggsave(
    file.path(OUT_FIGS, paste0("FIG_D_AGE_SPECIFIC_DIFF_", i, "_", TABLE_PICK, "_", ROI_PICK, ".png")),
    p_diff, width = 7.0, height = 5.3, dpi = DPI_FIG, bg = "white"
  )
}

################################################################################
# EXTRA BLOCK — FINAL FIGURES FOR AGE × cPES RESULTS
#
# Includes:
#   FIG_E1 = 2D continuous prediction surface
#   FIG_E2 = 2D surface + observed points
#   FIG_E3 = 3D static perspective
#   FIG_E4 = slope heatmap
#   FIG_E5 = 3D wireframe static
#   FIG_F1 = scatterplot with all participants
#   FIG_F2 = scatterplot with low/high cPES phenotype only
################################################################################

build_lmm_surface_grid <- function(model_obj,
                                   age_seq,
                                   cpes_seq,
                                   sex_value = "Male",
                                   site_value = NULL,
                                   site_levels = NULL) {
  
  nd <- expand.grid(
    age_wave = age_seq,
    cPES_z = cpes_seq,
    gender = factor(sex_value, levels = c("Male", "Female")),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  model_terms <- attr(terms(model_obj), "term.labels")
  if ("site_original" %in% model_terms) {
    nd$site_original <- factor(site_value, levels = site_levels)
  }
  
  X <- model.matrix(delete.response(terms(model_obj)), nd)
  beta <- lme4::fixef(model_obj)
  V <- as.matrix(vcov(model_obj))
  
  nd$predicted <- as.numeric(X %*% beta)
  nd$se_fit <- sqrt(diag(X %*% V %*% t(X)))
  nd$conf.low <- nd$predicted - 1.96 * nd$se_fit
  nd$conf.high <- nd$predicted + 1.96 * nd$se_fit
  
  nd
}

matrix_from_grid <- function(grid_df, value_col = "predicted") {
  x_vals <- sort(unique(grid_df$age_wave))
  y_vals <- sort(unique(grid_df$cPES_z))
  
  z_mat <- outer(
    x_vals, y_vals,
    Vectorize(function(x, y) {
      v <- grid_df[[value_col]][grid_df$age_wave == x & grid_df$cPES_z == y]
      if (length(v) == 0) NA_real_ else as.numeric(v[1])
    })
  )
  
  storage.mode(z_mat) <- "double"
  z_mat
}

compute_conditional_age_slope_full <- function(model_obj, cpes_values, sex = "Male"){
  tb <- tidy(model_obj, effects = "fixed", conf.int = FALSE)
  
  get_coef <- function(term){
    val <- tb$estimate[tb$term == term][1]
    if (length(val) == 0 || is.na(val)) 0 else val
  }
  
  vc <- as.matrix(vcov(model_obj))
  rn <- rownames(vc)
  
  b_age      <- get_coef("age_wave")
  b_axe      <- get_coef("age_wave:cPES_z")
  b_age_fem  <- get_coef("age_wave:genderFemale")
  b_axesex   <- get_coef("age_wave:cPES_z:genderFemale")
  
  female_ind <- ifelse(sex == "Female", 1, 0)
  
  out <- lapply(cpes_values, function(e){
    slope <- b_age + b_axe * e + b_age_fem * female_ind + b_axesex * e * female_ind
    
    coef_vec <- rep(0, length(rn))
    names(coef_vec) <- rn
    if ("age_wave" %in% rn) coef_vec["age_wave"] <- 1
    if ("age_wave:cPES_z" %in% rn) coef_vec["age_wave:cPES_z"] <- e
    if (female_ind == 1 && "age_wave:genderFemale" %in% rn) coef_vec["age_wave:genderFemale"] <- 1
    if (female_ind == 1 && "age_wave:cPES_z:genderFemale" %in% rn) coef_vec["age_wave:cPES_z:genderFemale"] <- e
    
    se <- sqrt(as.numeric(t(coef_vec) %*% vc %*% coef_vec))
    
    tibble(
      cPES_z = e,
      sex = sex,
      slope_age = slope,
      conf.low = slope - 1.96 * se,
      conf.high = slope + 1.96 * se
    )
  })
  
  bind_rows(out)
}

save_surface_persp_png <- function(grid_df,
                                   file_png,
                                   title_txt,
                                   zlab_txt = "Predicted thickness") {
  
  x_vals <- sort(unique(grid_df$age_wave))
  y_vals <- sort(unique(grid_df$cPES_z))
  z_mat  <- matrix_from_grid(grid_df, "predicted")
  
  z_mid <- (
    z_mat[-1, -1] +
      z_mat[-nrow(z_mat), -1] +
      z_mat[-1, -ncol(z_mat)] +
      z_mat[-nrow(z_mat), -ncol(z_mat)]
  ) / 4
  
  pal <- colorRampPalette(c(COL_LOW, "grey95", COL_HIGH))(100)
  z_scaled <- (z_mid - min(z_mid, na.rm = TRUE)) /
    (max(z_mid, na.rm = TRUE) - min(z_mid, na.rm = TRUE) + 1e-12)
  facet_cols <- pal[pmax(1, pmin(100, round(z_scaled * 99) + 1))]
  
  png(file_png, width = 2600, height = 2100, res = 320, bg = "white")
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  
  par(mar = c(4.2, 4.2, 3.0, 2.0))
  
  persp(
    x = x_vals,
    y = y_vals,
    z = z_mat,
    theta = 42,
    phi = 24,
    expand = 0.9,
    col = facet_cols,
    border = NA,
    shade = 0.35,
    ticktype = "detailed",
    xlab = "Age (years)",
    ylab = "cPES z",
    zlab = zlab_txt,
    main = title_txt,
    cex.lab = 1.0,
    cex.axis = 0.9
  )
}

save_surface_wire_png <- function(grid_df,
                                  file_png,
                                  title_txt,
                                  zlab_txt = "Predicted thickness") {
  
  x_vals <- sort(unique(grid_df$age_wave))
  y_vals <- sort(unique(grid_df$cPES_z))
  z_mat  <- matrix_from_grid(grid_df, "predicted")
  
  png(file_png, width = 2600, height = 2100, res = 320, bg = "white")
  op <- par(no.readonly = TRUE)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)
  
  par(mar = c(4.2, 4.2, 3.0, 2.0))
  
  pm <- persp(
    x = x_vals,
    y = y_vals,
    z = z_mat,
    theta = 34,
    phi = 20,
    expand = 0.85,
    col = "white",
    border = "grey35",
    shade = 0,
    ticktype = "detailed",
    xlab = "Age (years)",
    ylab = "cPES z",
    zlab = zlab_txt,
    main = title_txt,
    cex.lab = 1.0,
    cex.axis = 0.9
  )
  
  for (ii in seq_len(nrow(z_mat))) {
    xy <- trans3d(
      x = rep(x_vals[ii], length(y_vals)),
      y = y_vals,
      z = z_mat[ii, ],
      pmat = pm
    )
    lines(xy, col = "grey20", lwd = 0.7)
  }
  
  for (jj in seq_len(ncol(z_mat))) {
    xy <- trans3d(
      x = x_vals,
      y = rep(y_vals[jj], length(x_vals)),
      z = z_mat[, jj],
      pmat = pm
    )
    lines(xy, col = "grey20", lwd = 0.7)
  }
}

if (nrow(SHORT_AXE) > 0) {
  
  for (i in seq_len(nrow(SHORT_AXE))) {
    
    TABLE_PICK <- SHORT_AXE$table[i]
    ROI_PICK   <- SHORT_AXE$roi[i]
    ROI_LABEL  <- SHORT_AXE$roi_label[i]
    
    wide_pick <- read_csv(
      file.path(OUT_QC, paste0("00_", TABLE_PICK, "_wide_merged.csv")),
      show_col_types = FALSE
    ) %>% clean_names()
    
    roi_col <- resolve_roi_col(wide_pick, ROI_PICK)
    cpes_col <- get_cpes_col_after_clean(wide_pick)
    if (is.na(cpes_col)) stop("cPES_z / cpes_z column not found after clean_names() in: ", TABLE_PICK)
    
    dat_pick <- wide_pick %>%
      transmute(
        ident = as.character(ident),
        age_wave = as.numeric(age_wave),
        cPES_z = as.numeric(.data[[cpes_col]]),
        gender = factor(gender, levels = c("Male", "Female")),
        site_original = factor(site_original),
        value = suppressWarnings(as.numeric(.data[[roi_col]]))
      ) %>%
      filter(!is.na(value), !is.na(age_wave), !is.na(cPES_z), !is.na(gender))
    
    use_site_pick <- USE_SITE && ("site_original" %in% names(dat_pick)) &&
      (nlevels(dat_pick$site_original) >= 2)
    
    ref_site <- if (use_site_pick) pick_reference_site(dat_pick$site_original) else NA_character_
    
    f_plot <- as.formula(paste0(
      "value ~ age_wave * cPES_z * gender",
      if (use_site_pick) " + site_original" else "",
      " + (1 | ident)"
    ))
    
    fit_plot <- lmer(
      f_plot,
      data = dat_pick,
      REML = FALSE,
      control = lmerControl(calc.derivs = FALSE)
    )
    
    age_seq <- seq(min(dat_pick$age_wave, na.rm = TRUE),
                   max(dat_pick$age_wave, na.rm = TRUE),
                   length.out = 120)
    
    cpes_seq <- seq(quantile(dat_pick$cPES_z, 0.05, na.rm = TRUE),
                   quantile(dat_pick$cPES_z, 0.95, na.rm = TRUE),
                   length.out = 120)
    
    surf_df <- build_lmm_surface_grid(
      model_obj = fit_plot,
      age_seq = age_seq,
      cpes_seq = cpes_seq,
      sex_value = "Male",
      site_value = ref_site,
      site_levels = if (use_site_pick) levels(dat_pick$site_original) else NULL
    ) %>%
      mutate(
        table = TABLE_PICK,
        roi = ROI_PICK,
        roi_label = ROI_LABEL,
        sex_fixed = "Male",
        site_fixed = ifelse(use_site_pick, ref_site, NA_character_)
      )
    
    write_csv(
      surf_df,
      file.path(OUT_TAB, paste0("15_SURFACE_GRID_CONTINUOUS_", i, "_", TABLE_PICK, "_", ROI_PICK, ".csv"))
    )
    
    # FIG_E1 — continuous 2D surface
    p_surface2d <- ggplot(surf_df, aes(x = age_wave, y = cPES_z)) +
      geom_raster(aes(fill = predicted), interpolate = TRUE) +
      geom_contour(aes(z = predicted), color = "white", alpha = 0.55, linewidth = 0.30) +
      theme_paper(base_size = 12) +
      scale_fill_gradient2(
        low = COL_LOW,
        mid = "grey95",
        high = COL_HIGH,
        midpoint = median(surf_df$predicted, na.rm = TRUE),
        name = "Predicted\nthickness"
      ) +
      labs(
        x = "Age (years)",
        y = "Exposome Risk Score (cPES z)",
        title = paste0("Continuous age × cPES effect: ", ROI_LABEL),
        subtitle = "Population-level predictions from the primary LMM"
      ) +
      theme(
        legend.position = "right",
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(
      file.path(OUT_FIGS, paste0("FIG_E1_SURFACE2D_CONTINUOUS_", i, "_", TABLE_PICK, "_", ROI_PICK, ".png")),
      p_surface2d, width = 7.4, height = 5.9, dpi = DPI_FIG, bg = "white"
    )
    
    # FIG_E2 — continuous 2D surface + observed points
    dat_points <- dat_pick %>%
      filter(gender == "Male") %>%
      mutate(
        age_wave_j = jitter(age_wave, amount = 0.05),
        cPES_z_j = jitter(cPES_z, amount = 0.03)
      )
    
    p_surface_points <- ggplot() +
      geom_raster(
        data = surf_df,
        aes(x = age_wave, y = cPES_z, fill = predicted),
        interpolate = TRUE
      ) +
      geom_contour(
        data = surf_df,
        aes(x = age_wave, y = cPES_z, z = predicted),
        color = "white",
        alpha = 0.55,
        linewidth = 0.28
      ) +
      geom_point(
        data = dat_points,
        aes(x = age_wave_j, y = cPES_z_j),
        inherit.aes = FALSE,
        shape = 21,
        fill = "black",
        color = "black",
        alpha = 0.35,
        size = 1.2,
        stroke = 0.1
      ) +
      theme_paper(base_size = 12) +
      scale_fill_gradient2(
        low = COL_LOW,
        mid = "grey95",
        high = COL_HIGH,
        midpoint = median(surf_df$predicted, na.rm = TRUE),
        name = "Predicted\nthickness"
      ) +
      labs(
        x = "Age (years)",
        y = "Exposome Risk Score (cPES z)",
        title = paste0("Observed data over continuous age × cPES surface: ", ROI_LABEL),
        subtitle = "Black points represent observed male/reference-level data"
      ) +
      theme(
        legend.position = "right",
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(
      file.path(OUT_FIGS, paste0("FIG_E2_SURFACE2D_POINTS_", i, "_", TABLE_PICK, "_", ROI_PICK, ".png")),
      p_surface_points, width = 7.4, height = 5.9, dpi = DPI_FIG, bg = "white"
    )
    
    # FIG_E3 — 3D static perspective
    save_surface_persp_png(
      grid_df = surf_df,
      file_png = file.path(OUT_FIGS, paste0("FIG_E3_SURFACE3D_CONTINUOUS_", i, "_", TABLE_PICK, "_", ROI_PICK, ".png")),
      title_txt = paste0("Continuous age × cPES surface: ", ROI_LABEL),
      zlab_txt = "Predicted thickness"
    )
    
    # FIG_E4 — slope heatmap
    slope_df <- compute_conditional_age_slope_full(
      model_obj = fit_plot,
      cpes_values = cpes_seq,
      sex = "Male"
    ) %>%
      mutate(
        age_wave = mean(age_seq),
        roi_label = ROI_LABEL
      )
    
    slope_heat_df <- expand.grid(
      age_wave = age_seq,
      cPES_z = cpes_seq,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    ) %>%
      left_join(
        slope_df %>% select(cPES_z, slope_age, conf.low, conf.high),
        by = "cPES_z"
      ) %>%
      mutate(roi_label = ROI_LABEL)
    
    write_csv(
      slope_heat_df,
      file.path(OUT_TAB, paste0("16_SLOPE_HEATMAP_DATA_", i, "_", TABLE_PICK, "_", ROI_PICK, ".csv"))
    )
    
    max_abs_slope <- max(abs(slope_heat_df$slope_age), na.rm = TRUE)
    if (!is.finite(max_abs_slope) || is.na(max_abs_slope)) max_abs_slope <- 0.01
    
    p_slope_heat <- ggplot(slope_heat_df, aes(x = age_wave, y = cPES_z)) +
      geom_raster(aes(fill = slope_age), interpolate = TRUE) +
      geom_contour(aes(z = slope_age), color = "white", alpha = 0.55, linewidth = 0.28) +
      theme_paper(base_size = 12) +
      scale_fill_gradient2(
        low = COL_HIGH,
        mid = "grey95",
        high = COL_LOW,
        midpoint = 0,
        limits = c(-max_abs_slope, max_abs_slope),
        name = "Age slope\n(mm/year)"
      ) +
      labs(
        x = "Age (years)",
        y = "Exposome Risk Score (cPES z)",
        title = paste0("Conditional age slope across cPES: ", ROI_LABEL),
        subtitle = "This plot directly represents the age × cPES interaction term"
      ) +
      theme(
        legend.position = "right",
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(
      file.path(OUT_FIGS, paste0("FIG_E4_SLOPE_HEATMAP_", i, "_", TABLE_PICK, "_", ROI_PICK, ".png")),
      p_slope_heat, width = 7.4, height = 5.9, dpi = DPI_FIG, bg = "white"
    )
    
    # FIG_E5 — 3D wireframe static
    save_surface_wire_png(
      grid_df = surf_df,
      file_png = file.path(OUT_FIGS, paste0("FIG_E5_SURFACE3D_WIREFRAME_", i, "_", TABLE_PICK, "_", ROI_PICK, ".png")),
      title_txt = paste0("Wireframe age × cPES surface: ", ROI_LABEL),
      zlab_txt = "Predicted thickness"
    )
    
    # FIG_F1 — scatterplot with all participants
    p_scatter_all <- ggplot(dat_pick, aes(x = age_wave, y = value, color = cPES_z)) +
      geom_point(alpha = 0.45, size = 1.5) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black", linewidth = 0.8) +
      theme_paper(base_size = 12) +
      scale_color_gradient2(
        low = COL_LOW,
        mid = "grey70",
        high = COL_HIGH,
        midpoint = 0,
        name = "cPES z"
      ) +
      labs(
        x = "Age (years)",
        y = "Cortical thickness",
        title = paste0(ROI_LABEL),
        subtitle = "Scatterplot colored by continuous cPES"
      ) +
      theme(
        legend.position = "right",
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(
      file.path(OUT_FIGS, paste0("FIG_F1_SCATTER_ALL_", i, "_", TABLE_PICK, "_", ROI_PICK, ".png")),
      p_scatter_all, width = 7.2, height = 5.6, dpi = DPI_FIG, bg = "white"
    )
    
    # FIG_F2 — scatterplot with low/high phenotype only
    dat_lowhigh <- dat_pick %>%
      mutate(
        cPES_group = case_when(
          cPES_z <= -1 ~ "Low cPES",
          cPES_z >=  1 ~ "High cPES",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(cPES_group)) %>%
      mutate(
        cPES_group = factor(cPES_group, levels = c("Low cPES", "High cPES"))
      )
    
    write_csv(
      dat_lowhigh,
      file.path(OUT_TAB, paste0("17_LOW_HIGH_PHENOTYPE_DATA_", i, "_", TABLE_PICK, "_", ROI_PICK, ".csv"))
    )
    
    n_low_sc <- dat_lowhigh %>% distinct(ident, cPES_group) %>% filter(cPES_group == "Low cPES") %>% nrow()
    n_high_sc <- dat_lowhigh %>% distinct(ident, cPES_group) %>% filter(cPES_group == "High cPES") %>% nrow()
    
    dat_lowhigh <- dat_lowhigh %>%
      mutate(
        cPES_group_n = factor(
          cPES_group,
          levels = c("Low cPES", "High cPES"),
          labels = c(
            paste0("Low cPES (n=", n_low_sc, ")"),
            paste0("High cPES (n=", n_high_sc, ")")
          )
        )
      )
    
    p_scatter_lowhigh <- ggplot(dat_lowhigh, aes(x = age_wave, y = value, color = cPES_group_n)) +
      geom_point(alpha = 0.55, size = 1.7) +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.9) +
      theme_paper(base_size = 12) +
      scale_color_manual(
        values = setNames(
          c(COL_LOW, COL_HIGH),
          c(
            paste0("Low cPES (n=", n_low_sc, ")"),
            paste0("High cPES (n=", n_high_sc, ")")
          )
        ),
        name = NULL
      ) +
      labs(
        x = "Age (years)",
        y = "Cortical thickness",
        title = paste0("Observed low/high cPES phenotype: ", ROI_LABEL),
        subtitle = "Phenotype defined as cPES ≤ -1 SD vs cPES ≥ +1 SD"
      ) +
      theme(
        legend.position = "bottom",
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(
      file.path(OUT_FIGS, paste0("FIG_F2_SCATTER_LOWHIGH_", i, "_", TABLE_PICK, "_", ROI_PICK, ".png")),
      p_scatter_lowhigh, width = 7.2, height = 5.6, dpi = DPI_FIG, bg = "white"
    )
  }
  
} else {
  msg("  No age × cPES cortical ROIs -> skipping final continuous/scatter plots.")
}


