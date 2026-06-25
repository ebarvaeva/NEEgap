# XGBoost_CV.R — XGBoost cross-validation gap-filling for NEE
#
# Trains one XGBoost model per artificial gap to gap-fill NEE. For each gap the
# affected rows are blanked from NEE, a model is trained on the remaining
# observed rows (with a stratified train/validation split), and the blanked rows
# are predicted back — an out-of-sample test of gap-filling skill across four
# gap-size classes (VL = very large, L = large, M = medium, S = small).
#
# Inputs (provided by the run script, assumed to already exist):
#   df          — JCi_cv.rds: observed-NEE rows, pre-built by Scripts 08-09 with
#                 gap-flag columns (VL1…, L1…, M1…, S1…), masked-NEE columns
#                 (NEE_VL1…), and PI columns (PI_VL1…). Needs timestamp + night
#                 columns, used to stratify the validation split.
#   predictors  — character vector of base predictor column names
#   RESULTS_DIR — output directory (already created)
#   PI_ENABLED  — logical; if TRUE, the gap-specific PI_{label} column is added
#                 to the predictors for that gap
#
# Outputs (written to RESULTS_DIR):
#   df_cv_all_predictions.rds              — gap-filled NEE predictions
#   boosting_evaluation_curves/
#     XGB_NEE_{SIZE}_{LABEL}_rmse.png      — train/validation RMSE per round


# fix the seed and pin BLAS to one thread so runs are reproducible
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

# stop early with a clear message if any required package is missing
needed <- c("dplyr", "tibble", "tidyr", "ggplot2", "lubridate", "xgboost")
miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))

# month -> season label (Nov-Jan = Winter, Feb-Apr = Spring, …)
timestamp_to_season <- function(x) {
  m <- as.integer(format(x, "%m"))
  dplyr::case_when(m %in% c(11,12,1) ~ "Winter", m %in% c(2,3,4) ~ "Spring",
                   m %in% c(5,6,7)   ~ "Summer", m %in% c(8,9,10) ~ "Autumn",
                   TRUE ~ NA_character_)
}

# stratified train/validation split, balanced across season|month and day/night
stratified_train_val_split <- function(flux_data, obs_rows, val_frac = 0.10, seed = 42L) {
  n_total <- length(obs_rows)
  if (!n_total) return(list(train_rows = integer(0), val_rows = integer(0)))
  n_val <- max(1L, round(val_frac * n_total))
  
  # parse timestamps to season + month, and read the day/night flag
  ts <- flux_data$timestamp[obs_rows]
  if (inherits(ts, "Date"))    ts <- as.POSIXct(ts)
  if (!inherits(ts, "POSIXt")) ts <- suppressWarnings(lubridate::parse_date_time(ts, orders = c("Ymd HMS","Ymd HM","Ymd")))
  seas <- timestamp_to_season(ts)
  mon  <- as.integer(format(ts, "%m"))
  ngt  <- suppressWarnings(as.integer(flux_data$night[obs_rows]))
  
  # keep only rows with a usable season/month/flag
  ok <- !(is.na(seas) | is.na(mon) | is.na(ngt))
  idx_ok <- obs_rows[ok]; sk <- seas[ok]; mk <- mon[ok]; nk <- ngt[ok]
  
  set.seed(seed)
  n0 <- floor(n_val / 2L); n1 <- n_val - n0   # split validation between day and night
  pick <- integer(0)
  
  # round-robin: draw `target` rows as evenly as possible across the strata pools
  rr <- function(pools, target) {
    if (!length(pools) || target <= 0) return(integer(0))
    szs <- vapply(pools, length, integer(1)); q <- integer(length(pools)); left <- target; ord <- order(names(pools))
    while (left > 0) {
      prog <- FALSE
      for (j in ord) { if (left <= 0) break; if (q[j] < szs[j]) { q[j] <- q[j] + 1L; left <- left - 1L; prog <- TRUE } }
      if (!prog) break
    }
    unlist(Map(function(ids, k) if (k > 0L) sample(ids, k) else integer(0), pools, as.list(q)), use.names = FALSE)
  }
  
  # draw day rows and night rows from their season|month strata
  if (any(ok) && any(nk == 0L)) pick <- c(pick, rr(split(idx_ok[nk == 0L], paste(sk[nk == 0L], mk[nk == 0L], sep = "|")), n0))
  if (any(ok) && any(nk == 1L)) pick <- c(pick, rr(split(idx_ok[nk == 1L], paste(sk[nk == 1L], mk[nk == 1L], sep = "|")), n1))
  
  # top up to n_val if the strata could not supply enough
  sf <- n_val - length(pick); if (sf > 0) { rem <- setdiff(idx_ok, pick);   if (length(rem)) pick <- c(pick, sample(rem, min(sf, length(rem)))) }
  sf <- n_val - length(pick); if (sf > 0) { rem <- setdiff(obs_rows, pick); if (length(rem)) pick <- c(pick, sample(rem, min(sf, length(rem)))) }
  
  list(train_rows = sort(setdiff(obs_rows, pick)), val_rows = sort(unique(pick)))
}

# build a numeric matrix (doubles) from selected rows and feature columns
make_xgb_matrix <- function(df, row_idx, cols) {
  X <- as.matrix(df[row_idx, cols, drop = FALSE])
  for (j in seq_along(cols)) X[, j] <- suppressWarnings(as.numeric(X[, j]))
  storage.mode(X) <- "double"; colnames(X) <- cols
  X
}

# fixed XGBoost hyperparameters (single-threaded, RMSE objective)
xgb_params <- list(objective = "reg:squarederror", eval_metric = "rmse", nthread = 1L)

# save the train/validation RMSE boosting curve for one gap as a PNG
save_xgb_boosting_curves <- function(trained_xgb, gap_size_cat, gap_label) {
  elog <- trained_xgb$evaluation_log
  if (is.null(elog) || !nrow(elog)) return(invisible(NULL))
  ct <- grep("^train_.*rmse$", names(elog), value = TRUE)
  ce <- grep("^eval_.*rmse$",  names(elog), value = TRUE)
  if (length(ct) != 1 || length(ce) != 1) return(invisible(NULL))
  
  # long-format RMSE-per-round for the two curves
  d <- tibble::tibble(round = seq_len(nrow(elog)),
                      `Training set`   = as.numeric(elog[[ct]]),
                      `Validation set` = as.numeric(elog[[ce]])) %>%
    tidyr::pivot_longer(-round, names_to = "dataset", values_to = "rmse")
  
  out_dir <- file.path(RESULTS_DIR, "boosting_evaluation_curves")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  p <- ggplot2::ggplot(d, ggplot2::aes(round, rmse, linetype = dataset)) +
    ggplot2::geom_line() + ggplot2::geom_point(size = 0.5) +
    ggplot2::labs(title = paste0("XGBoost NEE ", gap_size_cat, " '", gap_label, "' — boosting curve"),
                  subtitle = "Each round adds ONE new tree. RMSE = cumulative ensemble error.",
                  x = "Boosting Round", y = expression(RMSE~(mu*mol~m^{-2}~s^{-1})), linetype = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white"), legend.position = "top")
  ggplot2::ggsave(file.path(out_dir, sprintf("XGB_NEE_%s_%s_rmse.png", gap_size_cat, gap_label)),
                  p, width = 7.5, height = 4.5, dpi = 200)
  invisible(NULL)
}

# train one XGBoost per gap of a size class and write predictions into the gap rows
run_xgb_for_gap_size <- function(flux_data, gap_size_cat, feature_cols,
                                 val_frac = 0.10, n_rounds = 100L, seed = 42L) {
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
    
    # stratified train/validation split of the observed rows
    spl <- stratified_train_val_split(flux_data, obs_rows, val_frac, seed)
    tr <- spl$train_rows; vr <- spl$val_rows
    if (length(tr) < 5 || !length(vr)) next
    
    # train/validation/test matrices
    dtrain <- xgboost::xgb.DMatrix(data = make_xgb_matrix(flux_data, tr, feats),       label = as.numeric(flux_data[[masked_nee]][tr]), missing = NA_real_)
    dval   <- xgboost::xgb.DMatrix(data = make_xgb_matrix(flux_data, vr, feats),       label = as.numeric(flux_data[[masked_nee]][vr]), missing = NA_real_)
    dtest  <- xgboost::xgb.DMatrix(data = make_xgb_matrix(flux_data, gap_rows, feats), missing = NA_real_)
    
    # train, save the boosting curve, then predict the blanked rows
    set.seed(seed)
    xgb_model <- xgboost::xgb.train(params = xgb_params, data = dtrain, nrounds = as.integer(n_rounds),
                                    watchlist = list(train = dtrain, eval = dval), verbose = 0)
    save_xgb_boosting_curves(xgb_model, gap_size_cat, gap_lbl)
    flux_data[[pred_col]][gap_rows] <- as.numeric(predict(xgb_model, dtest))
  }
  flux_data
}

# df is pre-built by Scripts 08-09
flux_data <- df

# run every gap-size class in turn: very-large, large, medium, small
for (cat in c("VL", "L", "M", "S"))
  flux_data <- run_xgb_for_gap_size(flux_data, cat, predictors, n_rounds = 100L)

# save the gap-filled NEE predictions
saveRDS(flux_data, file.path(RESULTS_DIR, "df_cv_all_predictions.rds"))
