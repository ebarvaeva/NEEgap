# =============================================================================
# MDS (Marginal Distribution Sampling) — Cross-Validation Gap-Filling for NEE
# =============================================================================
#
# LOGIC: MDS gap-fills NEE only. Reco and GPP are not derived here.
#
# INPUTS
#   df — a pre-built JC1_cv.rds / JC2_cv.rds data frame, sourced externally.
#        Already filtered to rows with finite NEE_orig and augmented with
#        gap flag columns (VL1…, L1…, M1…, S1…) and masked NEE columns.
#
# OUTPUTS  (written to RESULTS_DIR)
#   df_cv_all_predictions.rds
#       NEE_{S|M|L|VL}_mds_predicted
#   progress.log | run_info.txt
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
needed <- c("dplyr", "tibble", "lubridate", "purrr")
miss   <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))
if (!requireNamespace("REddyProc", quietly = TRUE))
  stop("Package 'REddyProc' is required.  Run: install.packages('REddyProc')")


# =============================================================================
# SECTION 3 — Progress logger
# =============================================================================
if (!exists("RESULTS_DIR", inherits = TRUE) || is.null(RESULTS_DIR))
  RESULTS_DIR <- file.path(tempdir(), "mds_cv_fallback")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
.run_log <- list()
log_msg <- function(...) {
  txt <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste(..., collapse = ""))
  message(txt); .run_log <<- append(.run_log, list(txt))
  try(silent = TRUE, {
    writeLines(unlist(.run_log), file.path(RESULTS_DIR, "progress.log"))
    saveRDS(.run_log, file.path(RESULTS_DIR, "progress_log.rds"))
  }); invisible(txt)
}
log_msg("MDS (NEE gap-fill only) started.  RESULTS_DIR = ", RESULTS_DIR)


# =============================================================================
# SECTION 4 — Data validation
# =============================================================================
req_cols <- c("timestamp", "NEE_orig")
miss_cols <- setdiff(req_cols, names(df))
if (length(miss_cols)) stop("Missing required columns: ", paste(miss_cols, collapse = ", "))

if (!"VPD" %in% names(df)) {
  log_msg("Column 'VPD' not found — creating as NA (MDS will still run).")
  df$VPD <- NA_real_
}

if (!inherits(df$timestamp, "POSIXt"))
  df$timestamp <- suppressWarnings(
    lubridate::parse_date_time(df$timestamp,
                               orders = c("Ymd HMS", "Ymd HM", "Ymd", "Y/m/d HMS", "Y/m/d", "d/m/Y HMS", "d/m/Y"),
                               tz = "UTC"))
stopifnot(!any(is.na(df$timestamp)))

n_dup <- sum(duplicated(df$timestamp))
if (n_dup > 0) {
  log_msg(n_dup, " duplicated timestamps removed.")
  df <- df %>% dplyr::arrange(timestamp) %>% dplyr::distinct(timestamp, .keep_all = TRUE)
}
df <- dplyr::arrange(df, timestamp)
log_msg("Input rows: ", nrow(df), "  (pre-filtered to finite NEE_orig by Script 08)")


# =============================================================================
# SECTION 5 — REddyProc interface helpers
# =============================================================================
prepare_reddyproc_input <- function(source_df, target_col,
                                    ppfd_col = "PPFD", temp_col = "Temp", vpd_col = "VPD") {
  ts_min <- min(source_df$timestamp, na.rm = TRUE)
  ts_max <- max(source_df$timestamp, na.rm = TRUE)
  tibble::tibble(DateTime = seq(as.POSIXct(ts_min, tz = "UTC"),
                                as.POSIXct(ts_max, tz = "UTC"), by = "30 min")) %>%
    dplyr::left_join(
      source_df %>% dplyr::select(timestamp,
                                  NEE  = dplyr::all_of(target_col),
                                  Rg   = dplyr::all_of(ppfd_col),
                                  Tair = dplyr::all_of(temp_col),
                                  VPD  = dplyr::all_of(vpd_col)),
      by = c("DateTime" = "timestamp")) %>%
    dplyr::mutate(Year = lubridate::year(DateTime),
                  DoY  = lubridate::yday(DateTime),
                  Hour = lubridate::hour(DateTime) + lubridate::minute(DateTime) / 60)
}

run_one_mds <- function(eddy_in) {
  proc <- REddyProc::sEddyProc$new(
    ID       = "SITE",
    Data     = dplyr::select(eddy_in, DateTime, Year, DoY, Hour, NEE, Rg, Tair, VPD),
    ColNames = c("NEE", "Rg", "Tair", "VPD"), ColPOSIXTime = "DateTime", DTS = 48L)
  proc$sMDSGapFill("NEE")
  res      <- proc$sExportResults()
  fill_col <- intersect(c("NEE_fall", "NEE_f"), names(res))
  if (!length(fill_col)) stop("REddyProc export missing 'NEE_f' or 'NEE_fall'.")
  tibble::tibble(DateTime = eddy_in$DateTime, mds_prediction = as.numeric(res[[fill_col[1]]]))
}


# =============================================================================
# SECTION 6 — MDS cross-validation
# =============================================================================
# Two-pass design per gap label:
#   Pass 1 (baseline): MDS on NEE_orig — fills any real NAs in the output.
#   Pass 2 (per-gap):  MDS on NEE_{label} (NA at gap rows) — CV predictions.
# Predictions written ONLY to gap rows.

run_mds_for_gap_size <- function(df_cv, gap_size_cat) {
  log_msg("MDS [", gap_size_cat, "]: starting.")
  gap_labels <- names(df_cv)[grepl(paste0("^", gap_size_cat, "\\d+$"), names(df_cv))]
  gap_labels <- gap_labels[order(as.integer(sub(gap_size_cat, "", gap_labels)))]
  if (!length(gap_labels)) {
    log_msg("MDS [", gap_size_cat, "]: no gap columns — skipping."); return(df_cv)
  }
  log_msg("MDS [", gap_size_cat, "]: ", length(gap_labels), " gap labels.")
  
  pred_col          <- paste0("NEE_", gap_size_cat, "_mds_predicted")
  df_cv[[pred_col]] <- NA_real_
  
  baseline_preds <- tryCatch({
    eddy_base <- prepare_reddyproc_input(df_cv, target_col = "NEE_orig")
    res_base  <- run_one_mds(eddy_base)
    res_base$mds_prediction[match(df_cv$timestamp, res_base$DateTime)]
  }, error = function(e) {
    log_msg("MDS baseline failed: ", conditionMessage(e)); rep(NA_real_, nrow(df_cv))
  })
  real_na_rows <- which(!is.finite(df_cv$NEE_orig))
  if (length(real_na_rows))
    df_cv[[pred_col]][real_na_rows] <- baseline_preds[real_na_rows]
  
  for (gap_lbl in gap_labels) {
    masked_nee <- paste0("NEE_", gap_lbl)
    if (!masked_nee %in% names(df_cv)) next
    gap_rows <- which(df_cv[[gap_lbl]] %in% c(TRUE, 1))
    if (!length(gap_rows)) next
    
    gap_preds <- tryCatch({
      eddy_in <- prepare_reddyproc_input(df_cv, target_col = masked_nee)
      res     <- run_one_mds(eddy_in)
      res$mds_prediction[match(df_cv$timestamp, res$DateTime)]
    }, error = function(e) {
      log_msg("MDS failed for ", gap_lbl, ": ", conditionMessage(e))
      rep(NA_real_, nrow(df_cv))
    })
    
    df_cv[[pred_col]][gap_rows] <- gap_preds[gap_rows]
    log_msg("MDS [", gap_size_cat, "] ", gap_lbl, ": filled ", length(gap_rows), " gap rows.")
  }
  
  log_msg("MDS [", gap_size_cat, "] complete.")
  df_cv
}


# =============================================================================
# SECTION 7 — Run MDS for all gap sizes
# =============================================================================
log_msg("=== Running MDS cross-validation — NEE ===")
for (cat in c("VL", "L", "M", "S"))
  df <- run_mds_for_gap_size(df, cat)
log_msg("All MDS gap sizes complete.")


# =============================================================================
# SECTION 8 — Save output
# =============================================================================
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(df, file.path(RESULTS_DIR, "df_cv_all_predictions.rds"))
log_msg("Saved: df_cv_all_predictions.rds")
writeLines(c(
  paste("Run finished    :", as.character(Sys.time())),
  paste("R version       :", R.version.string),
  paste("Source RDS      :", rds_name),
  paste("RESULTS_DIR     :", RESULTS_DIR),
  paste("n_rows          :", nrow(df)),
  paste("Target          :", "NEE only (MDS)")
), file.path(RESULTS_DIR, "run_info.txt"))
log_msg("Saved: run_info.txt — MDS script finished.")
# ================================ end MDS_CV.R ================================