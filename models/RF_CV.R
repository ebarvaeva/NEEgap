# RF_CV.R — Random Forest cross-validation gap-filling for NEE
#
# Trains one Random Forest ({ranger}) per artificial gap to gap-fill NEE.
# For each gap the affected rows are blanked from NEE, an RF is trained on the
# remaining observed rows, and the blanked rows are predicted back — an
# out-of-sample test of gap-filling skill across four gap-size classes
# (VL = very large, L = large, M = medium, S = small).
#
# Inputs (provided by run_model.R, assumed to already exist):
#   df          — JCi_cv.rds: observed-NEE rows, pre-built by Scripts 08-09 with
#                 gap-flag columns (VL1…, L1…, M1…, S1…), masked-NEE columns
#                 (NEE_VL1…), and PI columns (PI_VL1…)
#   predictors  — character vector of base predictor column names
#   RESULTS_DIR — output directory (already created)
#   PI_ENABLED  — logical; if TRUE, the gap-specific PI_{label} column is added
#                 to the predictors for that gap
#
# Outputs (written to RESULTS_DIR):
#   df_cv_all_predictions.rds                      — gap-filled NEE predictions
#   variable_importance/
#     rf_vi_{SIZE}_{LABEL}.rds                     — one gap's importance
#     rf_variable_importance_{SIZE}_all.{rds,csv}  — all gaps of a size class
#
# Feature scaling is not required: RF threshold splits are invariant to monotone
# transforms (Hastie et al. 2009, §15.2).


# fix the seed and pin BLAS to one thread so runs are reproducible
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

# stop early with a clear message if any required package is missing
needed <- c("dplyr", "tibble", "glue", "ranger")
miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))

# month -> season label (Nov-Jan = Winter, Feb-Apr = Spring, …)
timestamp_to_season <- function(x) {
  m <- as.integer(format(x, "%m"))
  dplyr::case_when(m %in% c(11,12,1) ~ "Winter", m %in% c(2,3,4) ~ "Spring",
                   m %in% c(5,6,7)   ~ "Summer", m %in% c(8,9,10) ~ "Autumn",
                   TRUE ~ NA_character_)
}

# train one RF per gap of a size class and write predictions into the gap rows
run_rf_for_gap_size <- function(flux_data,
                                gap_size_cat,
                                feature_cols,
                                random_seed     = 42L,
                                importance_type = "impurity") {
  
  # gap-flag columns for this size class (S1, S2, …), ordered numerically
  gap_labels <- names(flux_data)[grepl(paste0("^", gap_size_cat, "\\d+$"), names(flux_data))]
  gap_labels <- gap_labels[order(as.integer(sub(gap_size_cat, "", gap_labels)))]
  
  # output column that will hold this size class's predictions
  prediction_col              <- paste0("NEE_", gap_size_cat, "_rf_predicted")
  flux_data[[prediction_col]] <- NA_real_
  
  # subfolder this script creates to store variable importance
  vi_dir <- file.path(RESULTS_DIR, "variable_importance")
  dir.create(vi_dir, recursive = TRUE, showWarnings = FALSE)
  vi_list <- list()
  
  for (gap_lbl in gap_labels) {
    masked_nee <- paste0("NEE_", gap_lbl)   # NEE with this gap blanked out
    
    # base predictors, plus this gap's PI column when PI is enabled
    feats <- feature_cols
    if (PI_ENABLED) feats <- unique(c(feats, paste0("PI_", gap_lbl)))
    
    train_rows <- which(is.finite(flux_data[[masked_nee]]))    # observed rows -> training set
    gap_rows   <- which(flux_data[[gap_lbl]] %in% c(TRUE, 1))  # blanked rows  -> prediction targets
    if (length(train_rows) < 2 || !length(gap_rows)) next      # skip gaps too small to train/predict
    
    # training frame: observed NEE ~ predictors
    train_df <- cbind.data.frame(
      NEE = as.numeric(flux_data[[masked_nee]][train_rows]),
      as.data.frame(flux_data[train_rows, feats, drop = FALSE])
    )
    
    # fit the random forest (NEE ~ all predictors)
    rf_model <- ranger::ranger(
      NEE ~ ., data = train_df,
      seed = random_seed, num.threads = 1L,
      importance = importance_type,
      respect.unordered.factors = "order"   # treat the season factor correctly
    )
    
    # save this gap's variable importance (raw + relative, sorted)
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
      vi_list[[gap_lbl]] <- dplyr::mutate(vi_tidy, gap_size = gap_size_cat, gap_label = gap_lbl)
    }
    
    # fill the blanked rows with the RF's predictions
    flux_data[[prediction_col]][gap_rows] <- as.numeric(
      predict(rf_model, data = as.data.frame(flux_data[gap_rows, feats, drop = FALSE]))$predictions
    )
  }
  
  # combine all gaps' importance for this size class and save as rds + csv
  if (length(vi_list)) {
    vi_all <- dplyr::bind_rows(vi_list)
    saveRDS(vi_all, file.path(vi_dir, glue::glue("rf_variable_importance_{gap_size_cat}_all.rds")))
    utils::write.csv(vi_all,
                     file.path(vi_dir, glue::glue("rf_variable_importance_{gap_size_cat}_all.csv")),
                     row.names = FALSE)
  }
  
  flux_data
}

# df is pre-built by Scripts 08-09; add the categorical season predictor
flux_data <- df
flux_data$season <- factor(
  timestamp_to_season(flux_data$timestamp),
  levels = c("Winter", "Spring", "Summer", "Autumn")
)
predictors <- unique(c(predictors, "season"))

# run every gap-size class in turn: very-large, large, medium, small
for (cat in c("VL", "L", "M", "S"))
  flux_data <- run_rf_for_gap_size(flux_data, cat, predictors)

# save the gap-filled NEE predictions
saveRDS(flux_data, file.path(RESULTS_DIR, "df_cv_all_predictions.rds"))
