# =============================================================================
# RF_CV.R — Random Forest Cross-Validation Gap-Filling for NEE
# =============================================================================
#
# WHAT THIS SCRIPT DOES
#   Trains a Random Forest (R package {ranger}) to gap-fill NEE only.
#   Reco and GPP are derived separately by Script 10 (flux_partitioning.R),
#   which reads the NEE predictions saved here.
#
# PRE-REQUISITES
#   This script expects df to be the output of Script 08 + Script 09, i.e.
#   JCi_cv.rds loaded by run_model.R.  That file already contains:
#     - Only rows where NEE_orig is observed  (Script 08 filtered NA rows)
#     - Gap flag columns  VL1…VLk, L1…Lk, M1…Mk, S1…Sk  (Script 08)
#     - Masked-NEE columns  NEE_VL1…, NEE_L1…, NEE_M1…, NEE_S1…  (Script 08)
#     - PI columns  PI_VL1…, PI_L1…, PI_M1…, PI_S1…  (Script 09)
#   Sections for NA filtering, gap construction, and PI computation are
#   therefore no longer needed and have been removed.
#   The PI column for each gap label is selected automatically inside
#   run_rf_for_gap_size(): for gap S1 it uses PI_S1, for M3 it uses PI_M3.
#
# HOW TO USE
#   Source from run_model.R, which must define:
#     df          — loaded from data/data_prepared/JCi_cv.rds
#     predictors  — character vector of base predictor column names
#     RESULTS_DIR — output directory path
#     rds_name    — source file name (written to run_info.txt)
#     PI_ENABLED  — logical; passed from run_model.R user settings
#
# OUTPUTS  (written to RESULTS_DIR)
#   df_cv_all_predictions.rds        ← NEE predictions only; Reco/GPP added by Script 10
#   variable_importance/
#       rf_vi_{SIZE}_{LABEL}.rds  |  rf_variable_importance_{SIZE}_all.{rds,csv}
#   progress.log | progress_log.rds | run_info.txt
#
# NOTE ON FEATURE SCALING
#   RF threshold splits are invariant to monotone feature transformations.
#   Min-max scaling is NOT required.  Ref: Hastie et al. (2009) §15.2.
#
# =============================================================================


# =============================================================================
# SECTION 1 — Reproducibility
# =============================================================================

set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")


# =============================================================================
# SECTION 2 — Required packages
# =============================================================================

needed <- c("dplyr", "tibble", "tidyr", "ggplot2", "glue",
            "lubridate", "rlang", "purrr", "ranger")
miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))


# =============================================================================
# SECTION 3 — Progress logger
# =============================================================================

if (!exists("RESULTS_DIR", inherits = TRUE) || is.null(RESULTS_DIR))
  RESULTS_DIR <- file.path(tempdir(), "rf_cv_fallback")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

.run_log <- list()
log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
                paste(..., collapse = ""))
  message(txt); .run_log <<- append(.run_log, list(txt))
  try(silent = TRUE, {
    writeLines(unlist(.run_log), file.path(RESULTS_DIR, "progress.log"))
    saveRDS(.run_log, file.path(RESULTS_DIR, "progress_log.rds"))
  }); invisible(txt)
}
log_msg("RF (NEE gap-fill) started.  RESULTS_DIR = ", RESULTS_DIR)


# =============================================================================
# SECTION 4 — Helper functions
# =============================================================================

# 4a  Metrics
.vp <- function(o, p) is.finite(as.numeric(o)) & is.finite(as.numeric(p))
calc_mae  <- function(o, p) { ok<-.vp(o,p); if(!any(ok)) return(NA_real_); mean(abs(as.numeric(p)[ok]-as.numeric(o)[ok])) }
calc_rmse <- function(o, p) { ok<-.vp(o,p); if(!any(ok)) return(NA_real_); sqrt(mean((as.numeric(p)[ok]-as.numeric(o)[ok])^2)) }
calc_r2   <- function(o, p) { ok<-.vp(o,p); if(sum(ok)<2) return(NA_real_); ov<-as.numeric(o)[ok]; pv<-as.numeric(p)[ok]; den<-sum((ov-mean(ov))^2)*sum((pv-mean(pv))^2); if(den<=0) return(NA_real_); (sum((ov-mean(ov))*(pv-mean(pv)))^2)/den }

# 4b  Season labels
timestamp_to_season <- function(x) {
  m <- as.integer(format(x, "%m"))
  dplyr::case_when(m %in% c(11,12,1) ~ "Winter", m %in% c(2,3,4) ~ "Spring",
                   m %in% c(5,6,7)   ~ "Summer", m %in% c(8,9,10) ~ "Autumn",
                   TRUE ~ NA_character_)
}


# =============================================================================
# SECTION 5 — Phytomass Index flag
# =============================================================================
# PI_ENABLED is set by run_model.R and passed into this script via run_env.
# When TRUE, PI_{gap_label} columns (pre-computed by Script 09) are appended
# to the feature set for each gap: PI_S1 for gap S1, PI_M3 for gap M3, etc.
# When FALSE, no PI column is used and the base predictor set is unchanged.

if (!exists("PI_ENABLED")) PI_ENABLED <- FALSE   # safe fallback if sourced directly
log_msg("Phytomass Index: ", if (PI_ENABLED) "ENABLED" else "DISABLED")


# =============================================================================
# SECTION 6 — Random Forest cross-validation  (NEE only)
# =============================================================================
# For each gap label gl within gap_size_cat:
#   - Training rows = rows where NEE_{gl} is finite (i.e. outside the gap)
#   - If PI_ENABLED, PI_{gl} is appended to feats — the gap-specific PI
#     from Script 09, computed from outside-gl rows only
#   - One RF trained per gap; predictions written back to gap rows only

run_rf_for_gap_size <- function(flux_data,
                                gap_size_cat,
                                feature_cols,
                                random_seed     = 42L,
                                importance_type = "impurity") {
  log_msg("RF [", gap_size_cat, "]: starting.")
  if (!length(feature_cols)) stop("feature_cols is empty.")
  absent <- setdiff(feature_cols, names(flux_data))
  if (length(absent)) stop("Missing predictors: ", paste(absent, collapse=", "))
  
  gap_labels <- names(flux_data)[grepl(paste0("^", gap_size_cat, "\\d+$"), names(flux_data))]
  gap_labels <- gap_labels[order(as.integer(sub(gap_size_cat, "", gap_labels)))]
  if (!length(gap_labels)) {
    log_msg("RF [", gap_size_cat, "]: no gap columns — skipping."); return(flux_data)
  }
  log_msg("RF [", gap_size_cat, "]: ", length(gap_labels), " gaps to process.")
  
  prediction_col              <- paste0("NEE_", gap_size_cat, "_rf_predicted")
  flux_data[[prediction_col]] <- NA_real_
  
  vi_dir <- file.path(RESULTS_DIR, "variable_importance")
  dir.create(vi_dir, recursive = TRUE, showWarnings = FALSE)
  vi_list <- list()
  
  for (gap_lbl in gap_labels) {
    masked_nee <- paste0("NEE_", gap_lbl)
    if (!masked_nee %in% names(flux_data)) next
    
    feats <- feature_cols
    if (PI_ENABLED) {
      pi_col <- paste0("PI_", gap_lbl)   # PI_S1 for S1, PI_M3 for M3, etc.
      if (pi_col %in% names(flux_data)) feats <- unique(c(feats, pi_col))
    }
    
    train_rows <- which(is.finite(flux_data[[masked_nee]]))
    gap_rows   <- which(flux_data[[gap_lbl]] %in% c(TRUE, 1))
    log_msg("RF [", gap_size_cat, "] ", gap_lbl,
            ": n_train = ", length(train_rows), "  n_gap = ", length(gap_rows))
    if (length(train_rows) < 2 || !length(gap_rows)) next
    
    train_df <- cbind.data.frame(
      NEE = as.numeric(flux_data[[masked_nee]][train_rows]),
      as.data.frame(flux_data[train_rows, feats, drop = FALSE])
    )
    
    rf_model <- ranger::ranger(
      NEE ~ ., data = train_df,
      seed = random_seed, num.threads = 1L, importance = importance_type
    )
    
    if (!identical(importance_type, "none")) {
      vi_raw  <- rf_model$variable.importance
      vi_tidy <- tibble::tibble(
        predictor      = names(vi_raw),
        importance     = as.numeric(vi_raw),
        importance_rel = importance / sum(importance, na.rm = TRUE)
      ) %>% dplyr::arrange(dplyr::desc(importance))
      saveRDS(
        list(gap_size = gap_size_cat, gap_label = gap_lbl,
             n_train = length(train_rows), n_gap = length(gap_rows),
             feature_names = feats, importance_raw = vi_raw,
             importance_table = vi_tidy),
        file.path(vi_dir, glue::glue("rf_vi_{gap_size_cat}_{gap_lbl}.rds"))
      )
      vi_list[[gap_lbl]] <- dplyr::mutate(vi_tidy,
                                          gap_size  = gap_size_cat,
                                          gap_label = gap_lbl)
    }
    
    flux_data[[prediction_col]][gap_rows] <- as.numeric(
      predict(rf_model,
              data = as.data.frame(flux_data[gap_rows, feats, drop = FALSE]))$predictions
    )
  }
  
  if (length(vi_list)) {
    vi_all <- dplyr::bind_rows(vi_list)
    saveRDS(vi_all,
            file.path(vi_dir, glue::glue("rf_variable_importance_{gap_size_cat}_all.rds")))
    utils::write.csv(vi_all,
                     file.path(vi_dir, glue::glue("rf_variable_importance_{gap_size_cat}_all.csv")),
                     row.names = FALSE)
    log_msg("RF [", gap_size_cat, "]: variable importance saved.")
  }
  
  log_msg("RF [", gap_size_cat, "] complete — MAE = ",
          signif(calc_mae( flux_data$NEE_orig, flux_data[[prediction_col]]), 3),
          "  RMSE = ", signif(calc_rmse(flux_data$NEE_orig, flux_data[[prediction_col]]), 3),
          "  R² = ",   signif(calc_r2(  flux_data$NEE_orig, flux_data[[prediction_col]]), 3))
  flux_data
}


# =============================================================================
# SECTION 7 — Data preparation
# =============================================================================
# df arrives pre-filtered (no NA NEE_orig rows) and pre-built (all gap flag,
# masked-NEE, and PI columns already present from Scripts 08–09).

flux_data <- df
log_msg("Rows in pre-built CV data: ", nrow(flux_data))

# Report gap structure found in the pre-built data
for (pfx in c("VL","L","M","S")) {
  gcols <- names(flux_data)[grepl(paste0("^", pfx, "\\d+$"), names(flux_data))]
  if (length(gcols))
    log_msg("Gap columns found — ", pfx, ": ", length(gcols),
            " (", gcols[1], " … ", gcols[length(gcols)], ")")
}


# =============================================================================
# SECTION 8 — Run RF cross-validation  (NEE only)
# =============================================================================

log_msg("=== Running RF cross-validation — NEE ===")
for (cat in c("VL","L","M","S"))
  flux_data <- run_rf_for_gap_size(flux_data, cat, predictors)
log_msg("All RF NEE gap sizes complete.")


# =============================================================================
# SECTION 9 — Save outputs
# =============================================================================
# Reco and GPP are NOT derived here.  Run Script 10 (flux_partitioning.R)
# after this script to add Reco_{SIZE}_rf_predicted and GPP_{SIZE}_rf_predicted
# columns to df_cv_all_predictions.rds.

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(flux_data, file.path(RESULTS_DIR, "df_cv_all_predictions.rds"))
log_msg("Saved: df_cv_all_predictions.rds  (NEE predictions only — run Script 10 for Reco/GPP)")

writeLines(c(
  paste("Run finished    :", as.character(Sys.time())),
  paste("R version       :", R.version.string),
  paste("Platform        :", R.version$platform),
  paste("Source RDS      :", rds_name),
  paste("RESULTS_DIR     :", RESULTS_DIR),
  paste("n_rows          :", nrow(flux_data)),
  paste("Targets         :", "NEE only — Reco/GPP via Script 10"),
  paste("Predictors      :", paste(predictors, collapse = ", "))
), file.path(RESULTS_DIR, "run_info.txt"))
log_msg("Saved: run_info.txt — RF script finished.")

# ============================== end RF_CV.R ===================================