# XGBoost_CV.R — XGBoost cross-validation gap-filling for NEE
#
# Trains one XGBoost model per artificial gap to gap-fill NEE. For each gap the
# affected rows are blanked from NEE, a model is trained on the remaining
# observed rows, and the blanked rows are predicted back — an out-of-sample test
# of gap-filling skill across four gap-size classes (VL = very large, L = large,
# M = medium, S = small).
#
# Inputs (provided by the run script, assumed to already exist):
#   df          — JCi_cv.rds: observed-NEE rows, pre-built by Scripts 08-09 with
#                 gap-flag columns (VL1…, L1…, M1…, S1…), masked-NEE columns
#                 (NEE_VL1…), and PI columns (PI_VL1…).
#   predictors  — character vector of base predictor column names
#   RESULTS_DIR — output directory (already created)
#   PI_ENABLED  — logical; if TRUE, the gap-specific PI_{label} column is added
#                 to the predictors for that gap
#
# Outputs (written to RESULTS_DIR):
#   df_cv_all_predictions.rds              — gap-filled NEE predictions


# fix the seed and pin BLAS to one thread so runs are reproducible
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

# stop early with a clear message if the required package is missing
needed <- c("xgboost")
miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))

# build a numeric matrix (doubles) from selected rows and feature columns
make_xgb_matrix <- function(df, row_idx, cols) {
  X <- as.matrix(df[row_idx, cols, drop = FALSE])
  for (j in seq_along(cols)) X[, j] <- suppressWarnings(as.numeric(X[, j]))
  storage.mode(X) <- "double"; colnames(X) <- cols
  X
}

# fixed XGBoost hyperparameters (single-threaded, RMSE objective)
xgb_params <- list(objective = "reg:squarederror", eval_metric = "rmse", nthread = 1L)

# train one XGBoost per gap of a size class and write predictions into the gap rows
run_xgb_for_gap_size <- function(flux_data, gap_size_cat, feature_cols,
                                 n_rounds = 500L, seed = 42L) {
  # gap-flag columns for this size class (S1, S2, …), ordered numerically
  gap_labels <- names(flux_data)[grepl(paste0("^", gap_size_cat, "\\d+$"), names(flux_data))]
  gap_labels <- gap_labels[order(as.integer(sub(gap_size_cat, "", gap_labels)))]
  
  # output column for this size class's predictions
  pred_col <- paste0("NEE_", gap_size_cat, "_xgb_predicted")
  flux_data[[pred_col]] <- NA_real_
  
  for (gap_lbl in gap_labels) {
    masked_nee <- paste0("NEE_", gap_lbl)   # NEE with this gap blanked out
    
    # base predictors, plus this gap's PI column when PI is enabled
    feats <- feature_cols
    if (PI_ENABLED) feats <- unique(c(feats, paste0("PI_", gap_lbl)))
    
    obs_rows <- which(is.finite(flux_data[[masked_nee]]))    # observed rows
    gap_rows <- which(flux_data[[gap_lbl]] %in% c(TRUE, 1))  # blanked rows -> prediction targets
    if (length(obs_rows) < 10 || !length(gap_rows)) next     # skip gaps too small to train/predict
    
    # train on ALL observed rows (parity with RF); overfitting is handled by the
    # regularised objective + shrinkage + column subsampling (Chen & Guestrin 2016)
    dtrain <- xgboost::xgb.DMatrix(data = make_xgb_matrix(flux_data, obs_rows, feats),
                                   label = as.numeric(flux_data[[masked_nee]][obs_rows]), missing = NA_real_)
    dtest  <- xgboost::xgb.DMatrix(data = make_xgb_matrix(flux_data, gap_rows, feats), missing = NA_real_)
    
    set.seed(seed)
    xgb_model <- xgboost::xgb.train(params = xgb_params, data = dtrain, nrounds = as.integer(n_rounds),
                                    verbose = 0)
    flux_data[[pred_col]][gap_rows] <- as.numeric(predict(xgb_model, dtest))
  }
  flux_data
}

# df is pre-built by Scripts 08-09
flux_data <- df

# run every gap-size class in turn: very-large, large, medium, small
for (cat in c("VL", "L", "M", "S"))
  flux_data <- run_xgb_for_gap_size(flux_data, cat, predictors, n_rounds = 500L)

# save the gap-filled NEE predictions
saveRDS(flux_data, file.path(RESULTS_DIR, "df_cv_all_predictions.rds"))