# MDS_CV.R — Marginal Distribution Sampling cross-validation gap-filling for NEE
#
# Runs REddyProc MDS once per artificial gap to gap-fill NEE only (Reco and GPP
# are not derived here). For each gap the affected rows are blanked from NEE, MDS
# is run on the remaining series, and the blanked rows are predicted back — an
# out-of-sample test of gap-filling skill across four gap-size classes
# (VL = very large, L = large, M = medium, S = small).
#
# Inputs (provided by the run script, assumed to already exist):
#   df          — JCi_cv.rds: observed-NEE rows, pre-built by Scripts 08-09 with
#                 gap-flag columns (VL1…, L1…, M1…, S1…) and masked-NEE columns
#                 (NEE_VL1…)
#   RESULTS_DIR — output directory (already created)
#
# Outputs (written to RESULTS_DIR):
#   df_cv_all_predictions.rds — gap-filled NEE (NEE_{S|M|L|VL}_mds_predicted)


# fix the seed and pin BLAS to one thread so runs are reproducible
set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")

# stop early with a clear message if any required package is missing
needed <- c("dplyr", "tibble", "lubridate", "REddyProc")
miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))

# require timestamp + NEE_orig; add VPD as NA when absent (MDS still runs)
req_cols  <- c("timestamp", "NEE_orig")
miss_cols <- setdiff(req_cols, names(df))
if (length(miss_cols)) stop("Missing required columns: ", paste(miss_cols, collapse = ", "))
if (!"VPD" %in% names(df)) df$VPD <- NA_real_

# coerce timestamp to POSIXct, drop duplicates, sort
if (!inherits(df$timestamp, "POSIXt"))
  df$timestamp <- suppressWarnings(
    lubridate::parse_date_time(df$timestamp,
                               orders = c("Ymd HMS", "Ymd HM", "Ymd", "Y/m/d HMS", "Y/m/d", "d/m/Y HMS", "d/m/Y"),
                               tz = "UTC"))
stopifnot(!any(is.na(df$timestamp)))
df <- df %>% dplyr::arrange(timestamp) %>% dplyr::distinct(timestamp, .keep_all = TRUE)

# build a gap-free 30-min timeline and map df onto REddyProc's expected columns
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

# run one MDS gap-fill and return its NEE predictions aligned to DateTime
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

# per gap label: baseline MDS on NEE_orig fills real NAs, per-gap MDS fills the gap rows
run_mds_for_gap_size <- function(df_cv, gap_size_cat) {
  # gap-flag columns for this size class (S1, S2, …), ordered numerically
  gap_labels <- names(df_cv)[grepl(paste0("^", gap_size_cat, "\\d+$"), names(df_cv))]
  gap_labels <- gap_labels[order(as.integer(sub(gap_size_cat, "", gap_labels)))]
  if (!length(gap_labels)) return(df_cv)

  # output column that will hold this size class's predictions
  pred_col          <- paste0("NEE_", gap_size_cat, "_mds_predicted")
  df_cv[[pred_col]] <- NA_real_

  # baseline pass fills the genuinely missing NEE_orig rows
  baseline_preds <- tryCatch({
    eddy_base <- prepare_reddyproc_input(df_cv, target_col = "NEE_orig")
    res_base  <- run_one_mds(eddy_base)
    res_base$mds_prediction[match(df_cv$timestamp, res_base$DateTime)]
  }, error = function(e) rep(NA_real_, nrow(df_cv)))
  real_na_rows <- which(!is.finite(df_cv$NEE_orig))
  if (length(real_na_rows))
    df_cv[[pred_col]][real_na_rows] <- baseline_preds[real_na_rows]

  for (gap_lbl in gap_labels) {
    masked_nee <- paste0("NEE_", gap_lbl)   # NEE with this gap blanked out
    if (!masked_nee %in% names(df_cv)) next
    gap_rows <- which(df_cv[[gap_lbl]] %in% c(TRUE, 1))
    if (!length(gap_rows)) next

    # per-gap MDS predictions, kept only on the gap rows
    gap_preds <- tryCatch({
      eddy_in <- prepare_reddyproc_input(df_cv, target_col = masked_nee)
      res     <- run_one_mds(eddy_in)
      res$mds_prediction[match(df_cv$timestamp, res$DateTime)]
    }, error = function(e) rep(NA_real_, nrow(df_cv)))
    df_cv[[pred_col]][gap_rows] <- gap_preds[gap_rows]
  }
  df_cv
}

# run every gap-size class in turn: very-large, large, medium, small
for (cat in c("VL", "L", "M", "S"))
  df <- run_mds_for_gap_size(df, cat)

# save the gap-filled NEE predictions
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(df, file.path(RESULTS_DIR, "df_cv_all_predictions.rds"))
