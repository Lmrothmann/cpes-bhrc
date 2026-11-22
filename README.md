# Analysis scripts — cumulative psychosocial exposome score (cPES)

Essential, reproducible pipeline for "A cross-validated psychosocial exposome score predicts general
psychopathology and later psychiatric diagnosis in Brazilian children and adolescents" (BHRC).

## Data availability
Individual-level BHRC data are available through the cohort's data-access procedures (see cohort
profile, Salum … Pan, Int J Epidemiol 2025): these scripts document the full analytic procedure and
can be run by researchers who obtain data access. Personal file paths have been replaced by
placeholders (`<PROJECT_ROOT>`, `${cPES_OUT}`, `<BHRC_INPUT_DIR>`, `<COGNITION_DIR>`); set them to your
local tree before running. The fitted model coefficients (all 52 non-zero weights, eTable S2) are
provided so that the cPES itself is reproducible.

## Environment
R 4.4.1. Packages: tidyverse, haven, janitor, caret, glmnet, lme4, lmerTest, nlme, emmeans, broom,
broom.mixed, longCombat. Global random seed **241094** (set in every stochastic script).
Set the output root once: `Sys.setenv(cPES_OUT = "…/FINAL_WITH_SUBSTANCE")` (scripts default to a
relative `FINAL_WITH_SUBSTANCE`).

## Pipeline (run in order)
| # | script | produces |
|---|---|---|
| 00 | `00_data_curation.R` | curated exposure set (1,059 → 70) + cPES input matrix |
| 01 | `01_build_cPES.R` | nested-CV elastic net; **out-of-fold cPES**; full-sample coefficients (eTable S2); ExWAS (eFig1); selection stability (eFig2) |
| 02 | `02_longitudinal_cbcl_cognition.R` | linear mixed models of CBCL total and the 4 cognitive outcomes across waves (per-term estimates) |
| 03 | `03_cortical_thickness.R` | age×cPES×sex mixed models per Desikan-Killiany parcel (per-term estimates); male reference |
| 04a | `04a_cognition_adj_baselineCBCL.R` | cognitive models re-fitted with baseline CBCL as a covariate |
| 04b | `04b_cortical_adj_baselineCBCL.R` | cortical models re-fitted with baseline CBCL as a covariate |
| 05 | `05_site_combat_sensitivity.R` | site as fixed vs random effect and longitudinal ComBat, male reference (supersedes the buggy alphabetical-gender ComBat script) |
| 06 | `06_allmetrics_sweep.R` | age×cPES across area/volume/curvature/folding/subcortical (thickness-specificity check) |
| 07 | `07_domains.R` | descriptive 7-domain decomposition (in-sample) |
| 08 | `08_prospective_diagnosis.R` | logistic models of psychiatric diagnosis at follow-up, adjusted for baseline CBCL |
| 09 | `09_bootstrap_calibration.R` | bootstrap 95% CIs and calibration for out-of-fold performance |
| 11 | `11_sample_Ns_and_roi_table.R` | analytic sample sizes per wave (CBCL, each cognitive test, MRI) and the per-ROI table |
| 10 | `10_final_tables_corrected.R` | **canonical summary tables** (see below) — run LAST, after the model scripts |

Placeholders in every script (set before running): `<PROJECT_ROOT>`, `${cPES_OUT}`, `<BHRC_INPUT_DIR>`
(controlled BHRC .sav, not shareable), `<COGNITION_DIR>`.

## Important: multiple-comparison families are applied in script 10
The model scripts (02, 03) save correct per-term estimates (beta, SE, t, p). The **manuscript's
false-discovery-rate values and effect metric come from `10_final_tables_corrected.R`**, which applies:
- **cortex:** BH within hemisphere over the **34 Desikan-Killiany parcels only** (global measures eTIV,
  BrainSegVolNotVent and mean thickness are excluded) → the three right-hemisphere ROIs at q = 0.048;
- **cognition:** BH over the **four exposome-score coefficients only** (not stacked with intercepts /
  age terms) → three outcomes in the exploratory tier (q 0.064–0.080), none confirmatory;
- **developmental-years equivalent** with 95% CI.

Two-tier reporting convention: q < 0.05 = confirmatory, 0.05 ≤ q < 0.10 = exploratory (hypothesis-
generating). Exact q values are reported throughout.

## Coding notes (reproducibility)
- Sex enters cortical models as `factor(gender, levels = c("Male","Female"))` — **male is the reference
  category**, so `age_wave:cPES_z` is the male-referenced effect and the three-way `age×cPES×sex` term
  is the formal test of sex moderation. (Do not let `factor()` default to alphabetical order, which
  would silently make Female the reference and return the null effect.)
- Longitudinal CBCL/cognition models use assessment **wave** as the time metric with chronological age
  as a separate covariate; cortical models use continuous **age at scan** as the time metric.
- Excluded from this repo (not part of the final paper): CBCL-breadth / incremental-validity (circular),
  telomere analyses (companion paper), and a ComBat script with an alphabetical-gender bug — superseded
  by `05_site_combat_sensitivity.R`.
