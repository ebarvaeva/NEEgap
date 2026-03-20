# =============================================================================
# MDS_CV.R — Marginal Distribution Sampling (MDS) Cross-Validation Gap-Filling
#             NEE gap-fill with post-hoc Reco and GPP derivation
# =============================================================================
#
# PURPOSE
#   Fills artificial gaps in NEE using the Marginal Distribution Sampling
#   (MDS) algorithm implemented in the REddyProc package (Wutzler et al.,
#   2018, Biogeosciences).  MDS is a look-up table method: for each missing
#   half-hour, it searches a symmetric moving time window for half-hours with
#   similar radiation (Rg), temperature (Tair), and vapour pressure deficit
#   (VPD) conditions, and fills the gap with their mean NEE.  The window is
#   widened progressively (±5, ±10, ±20 days) until sufficient similar
#   conditions are found.
#
# MDS SERVES AS A BENCHMARK
#   MDS represents the standard process-independent gap-filling method used
#   in the eddy-covariance community (Reichstein et al., 2005; Papale et al.,
#   2006).  Including it alongside the machine-learning models allows direct
#   comparison with established practice.
#
# CROSS-VALIDATION ADAPTATION
#   REddyProc requires a complete, evenly-spaced half-hourly timeline without
#   NaN values in the forcing variables.  The cross-validation design therefore
#   operates differently from the ML models:
#   1. Gaps are constructed on the finite-NEE subset (flux_data).
#   2. Each gap label is processed by passing a complete annual data frame
#      to REddyProc with the gap rows' NEE set to NA.
#   3. REddyProc fills all missing values (genuine gaps + the current
#      artificial gap) using MDS.  Only the artificial gap rows are extracted
#      and written to the prediction columns.
#   This ensures REddyProc always operates on a realistic year-long dataset
#   rather than a stripped subset.
#
# INPUTS (passed from the calling run script)
#   df, RESULTS_DIR, rds_name
#   Note: predictors is not used by MDS.  The meteorological drivers
#   (Rg, Tair, VPD) are taken directly from the data frame columns.
#
# OUTPUTS  (written to RESULTS_DIR)
#   df_cv_all_predictions.rds
#     Original data frame augmented with:
#       NEE_{S|M|L|VL}_mds_predicted
#       Reco_{S|M|L|VL}_mds_predicted  (post-hoc partitioning)
#       GPP_{S|M|L|VL}_mds_predicted   (post-hoc partitioning)
#
#   progress.log | run_info.txt
#
# SECTIONS
#    1  Reproducibility            — random seed
#    2  Required packages          — includes REddyProc check
#    3  Progress logger
#    4  Data preparation           — full-year grid required by REddyProc
#    5  Utility functions          — gap partitioner, metrics, REddyProc adapter
#    6  REddyProc interface helpers — prepare_reddyproc_input(), run_one_mds()
#    7  Regrowth period detection  — for post-hoc Reco/GPP partitioning
#    8  Flux partitioning helpers  — same as RF_CV.R Section 7
#    9  MDS cross-validation       — run_mds_for_gap_size() definition
#   10  Construct artificial gaps
#   11  Run MDS cross-validation
#   12  Derive Reco and GPP
#   13  Save outputs
#
# REFERENCE
#   Wutzler et al. (2018) Basic and extensible post-processing of eddy
#   covariance flux data with EddyPro in combination with REddyProc.
#   Biogeosciences, 15, 5015-5030.
#
# =============================================================================


# =============================================================================
# SECTION 1 — Reproducibility
# =============================================================================
# Fixed seed and single-threaded execution.  See RF_CV.R Section 1.

set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")


# =============================================================================
# SECTION 2 — Required packages
# =============================================================================
# REddyProc is checked separately from the other packages because it has
# additional system-level dependencies and is not installed by default.

needed <- c("dplyr", "tibble", "tidyr", "lubridate", "purrr", "glue", "stringr")
miss   <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))
if (!requireNamespace("REddyProc", quietly = TRUE))
  stop("Package 'REddyProc' is required.  Run: install.packages('REddyProc')")


# =============================================================================
# SECTION 3 — Progress logger
# =============================================================================
# Timestamped logging to console and file.  See RF_CV.R Section 3.

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
log_msg("MDS (NEE gap-fill + post-hoc Reco/GPP) started.  RESULTS_DIR = ", RESULTS_DIR)


# =============================================================================
# SECTION 4 — Data preparation
# =============================================================================
# REddyProc requires a complete, evenly-spaced half-hourly timeline for the
# entire analysis period, with no missing time steps.  A full-year grid is
# constructed here and used as the base for all REddyProc calls.  Genuine
# measurement gaps (rows absent from the QC output) appear as NA in the
# forcing variables, which REddyProc can handle internally.
# Artificial cross-validation gaps are created on the finite-NEE subset
# (flux_data) and then mapped back onto the full-year grid for each gap label.
# This separation ensures that the gap-filling algorithm always sees a
# realistic year-long context, not a truncated subset.
# so REddyProc always receives a complete series.

req_cols <- c("timestamp", "NEE_orig", "PPFD", "Temp")
miss_cols <- setdiff(req_cols, names(df))
if (length(miss_cols)) stop("Missing required columns: ", paste(miss_cols, collapse = ", "))

if (!"VPD" %in% names(df)) {
  log_msg("Column 'VPD' not found — creating as NA (MDS will still run).")
  df$VPD <- NA_real_
}

if (!inherits(df$timestamp, "POSIXt"))
  df$timestamp <- suppressWarnings(
    lubridate::parse_date_time(df$timestamp,
      orders = c("Ymd HMS","Ymd HM","Ymd","Y/m/d HMS","Y/m/d","d/m/Y HMS","d/m/Y"),
      tz = "UTC"))
stopifnot(!any(is.na(df$timestamp)))

n_dup <- sum(duplicated(df$timestamp))
if (n_dup > 0) {
  log_msg(n_dup, " duplicated timestamps removed.")
  df <- df %>% dplyr::arrange(timestamp) %>% dplyr::distinct(timestamp, .keep_all = TRUE)
}
df <- dplyr::arrange(df, timestamp)

flux_data_finite <- df %>% dplyr::filter(is.finite(NEE_orig))
log_msg("Rows with finite NEE_orig (for gap construction): ", nrow(flux_data_finite))

RECO_GPP_ENABLED <- all(c("PPFD", "Temp") %in% names(flux_data_finite))
if (!RECO_GPP_ENABLED) log_msg("WARNING: PPFD or Temp missing — Reco/GPP will be skipped.")


# =============================================================================
# SECTION 5 — Utility functions
# =============================================================================
# Gap construction and metric computation.  create_cv_gaps_on_finite_subset()
# is a variant of create_cv_gaps() that builds gap matrices on the finite-NEE
# subset only, then returns index vectors that map back to the full data frame.
# This two-step approach is needed because REddyProc requires the full annual
# timeline, but the gap partitioning must operate only on rows with observed NEE.


# 5a  Block partitioner
partition_timeseries_into_blocks <- function(n, base_size, tol) {
  stopifnot(n > 0, base_size >= 1, tol >= 0)
  k <- max(1L, round(n / base_size)); sizes <- rep(base_size, k)
  delta <- n - sum(sizes)
  if (delta != 0) {
    m    <- min(k, max(1L, ceiling(abs(delta) / max(1, tol))))
    step <- sign(delta) * floor(abs(delta) / m)
    rem  <- sign(delta) * (abs(delta) - abs(step) * m)
    sizes[(k - m + 1):k] <- sizes[(k - m + 1):k] + step
    if (rem != 0) sizes[(k - abs(rem) + 1):k] <- sizes[(k - abs(rem) + 1):k] + sign(rem)
  }
  for (i in k:2) {
    dev <- sizes[i] - base_size
    if (abs(dev) > tol) {
      ex <- sign(dev) * (abs(dev) - tol)
      sizes[i] <- base_size + sign(dev) * tol; sizes[i - 1] <- sizes[i - 1] + ex
    }
  }
  sizes[sizes < 1] <- 1L; sizes[k] <- sizes[k] + (n - sum(sizes))
  starts <- c(1L, head(cumsum(sizes) + 1L, -1L)); ends <- cumsum(sizes)
  stopifnot(ends[length(ends)] == n)
  list(sizes = sizes, indices = purrr::map2(starts, ends, seq.int))
}

# 5b  Gap creator: builds flags on FINITE subset, joins back to full timeline
create_cv_gaps_on_finite_subset <- function(finite_data, full_data,
                                             base_size, tol, prefix) {
  n      <- nrow(finite_data)
  blocks <- partition_timeseries_into_blocks(n, as.integer(round(base_size)),
                                              as.integer(round(tol)))
  flag_names <- paste0(prefix, seq_along(blocks$indices))
  nee_names  <- paste0("NEE_", flag_names)
  flag_mat <- matrix(FALSE,    nrow=n, ncol=length(blocks$indices), dimnames=list(NULL,flag_names))
  nee_mat  <- matrix(NA_real_, nrow=n, ncol=length(blocks$indices), dimnames=list(NULL,nee_names))
  for (i in seq_along(blocks$indices)) {
    flag_mat[blocks$indices[[i]], i] <- TRUE
    nee_mat[, i] <- finite_data$NEE_orig
    nee_mat[blocks$indices[[i]], i] <- NA_real_
  }
  df_gaps   <- dplyr::bind_cols(finite_data,
                                 tibble::as_tibble(as.data.frame(flag_mat)),
                                 tibble::as_tibble(as.data.frame(nee_mat)))
  gap_meta  <- dplyr::select(df_gaps, timestamp,
                               dplyr::all_of(flag_names), dplyr::all_of(nee_names))
  full_joined <- dplyr::left_join(full_data, gap_meta, by = "timestamp") %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(flag_names),
                                ~ dplyr::coalesce(as.logical(.x), FALSE)))
  for (j in seq_along(flag_names))
    full_joined[[nee_names[j]]] <- ifelse(full_joined[[flag_names[j]]], NA_real_, full_joined$NEE_orig)
  list(data = dplyr::arrange(full_joined, timestamp), sizes = blocks$sizes)
}

# 5c  Metrics
.vp <- function(o,p) is.finite(as.numeric(o))&is.finite(as.numeric(p))
calc_mae  <- function(o,p){ok<-.vp(o,p); if(!any(ok)) return(NA_real_); mean(abs(as.numeric(p)[ok]-as.numeric(o)[ok]))}
calc_rmse <- function(o,p){ok<-.vp(o,p); if(!any(ok)) return(NA_real_); sqrt(mean((as.numeric(p)[ok]-as.numeric(o)[ok])^2))}
calc_r2   <- function(o,p){ok<-.vp(o,p); if(sum(ok)<2) return(NA_real_); ov<-as.numeric(o)[ok]; pv<-as.numeric(p)[ok]; den<-sum((ov-mean(ov))^2)*sum((pv-mean(pv))^2); if(den<=0) return(NA_real_); (sum((ov-mean(ov))*(pv-mean(pv)))^2)/den}


# =============================================================================
# SECTION 6 — REddyProc interface helpers
# =============================================================================
# Two functions that manage the REddyProc API:
#
# prepare_reddyproc_input(source_df, target_col, gap_rows, year)
#   Formats the full-year data frame for sEddyProc$sNew():
#   - Renames columns to the names expected by REddyProc (Rg, Tair, VPD, NEE)
#   - Sets the current artificial gap rows to NA in the NEE column
#   - Returns a data frame compatible with REddyProc's input format
#
# run_one_mds(eddy_in)
#   Creates a sEddyProc object, calls sMDSGapFill() for NEE, and extracts
#   the gap-filled NEE column.  MDS uses Rg, Tair, and VPD as look-up
#   variables with a maximum window of ±20 days.


prepare_reddyproc_input <- function(source_df, target_col,
                                     ppfd_col="PPFD", temp_col="Temp", vpd_col="VPD") {
  ts_min <- min(source_df$timestamp, na.rm=TRUE)
  ts_max <- max(source_df$timestamp, na.rm=TRUE)
  tibble::tibble(DateTime = seq(as.POSIXct(ts_min,tz="UTC"),
                                as.POSIXct(ts_max,tz="UTC"), by="30 min")) %>%
    dplyr::left_join(
      source_df %>% dplyr::select(timestamp,
                                   NEE  = dplyr::all_of(target_col),
                                   Rg   = dplyr::all_of(ppfd_col),
                                   Tair = dplyr::all_of(temp_col),
                                   VPD  = dplyr::all_of(vpd_col)),
      by = c("DateTime"="timestamp")) %>%
    dplyr::mutate(Year = lubridate::year(DateTime),
                  DoY  = lubridate::yday(DateTime),
                  Hour = lubridate::hour(DateTime) + lubridate::minute(DateTime)/60)
}

run_one_mds <- function(eddy_in) {
  proc <- REddyProc::sEddyProc$new(
    ID="SITE",
    Data=dplyr::select(eddy_in, DateTime,Year,DoY,Hour,NEE,Rg,Tair,VPD),
    ColNames=c("NEE","Rg","Tair","VPD"), ColPOSIXTime="DateTime", DTS=48L)
  proc$sMDSGapFill("NEE")
  res <- proc$sExportResults()
  fill_col <- intersect(c("NEE_fall","NEE_f"), names(res))
  if (!length(fill_col)) stop("REddyProc export missing 'NEE_f' or 'NEE_fall'.")
  tibble::tibble(DateTime=eddy_in$DateTime, mds_prediction=as.numeric(res[[fill_col[1]]]))
}


# =============================================================================
# SECTION 7 — Regrowth period detection
# =============================================================================
# Assigns regrowth_id per row from Grazing_days_since for use by the
# post-hoc flux partitioning in Section 8.  See RF_CV.R Section 6.

build_regrowth_periods <- function(days_since, threshold=1.0) {
  x <- suppressWarnings(as.numeric(days_since))
  is_zero <- is.finite(x)&(x==0); x_prev <- dplyr::lag(x)
  cumsum(ifelse(is_zero&(is.na(x_prev)|!is.finite(x_prev)|(x_prev>threshold)), 1L, 0L))
}


# =============================================================================
# SECTION 8 — Flux partitioning helpers
# =============================================================================
# Reco/GPP partitioning functions; identical to RF_CV.R Section 7.
# See that script for full algorithm documentation.
# Three functions shared by both ground-truth computation and gap prediction:
#
#  compute_reco_gpp_from_nee_orig()
#    Partitions NEE_orig into Reco_orig and GPP_orig per regrowth period.
#    These two columns serve as GROUND TRUTH when computing Reco/GPP metrics.
#    Called ONCE after regrowth periods are assigned.
#
#  derive_reco_gpp_from_filled_nee()
#    For each gap label Si, reconstructs a complete NEE series:
#      NEE_Si_with_predicted = NEE_Si  (observed outside gap)
#                             + NEE_{SIZE}_{model}_predicted  (at Si rows)
#    Then partitions this full series per regrowth period → Reco_Si, GPP_Si.
#    Writes only the Si rows into Reco_{SIZE}_{model}_predicted and
#    GPP_{SIZE}_{model}_predicted.
#    This ensures partitioning sees a complete, gap-free time series — the
#    same way standard post-hoc partitioning is applied in real pipelines.
#
# Night threshold : PPFD < 10 µmol m⁻² s⁻¹
# Day   threshold : PPFD > 10 µmol m⁻² s⁻¹
# Reco model      : Lloyd-Taylor Arrhenius  R10 estimated by OLS
# GPP  model      : Thornley non-rectangular hyperbola  fitted by BFGS

arrhenius_temp_scaling <- function(temp_celsius) {
  exp(309 * ((1 / (283.2 - 230)) - (1 / ((temp_celsius + 273.2) - 230))))
}

light_response_gpp <- function(params, ppfd_vec) {
  a <- params[1]; b <- params[2]; c <- params[3]
  num  <- a * ppfd_vec + b
  disc <- pmax(num^2 - 4 * c * (ppfd_vec * a * b), 0)
  (num - sqrt(disc)) / (2 * c)
}

gpp_ssr <- function(params, ppfd_vec, gpp_obs)
  sum((light_response_gpp(params, ppfd_vec) - gpp_obs)^2)

# Internal helper: fit R10 + GPP params and apply to a full period
.partition_one_period <- function(rows_in_period, nee_vec, ppfd_vec, temp_vec,
                                   is_night, is_day, gpp_init, ppfd_night_thr) {
  train_rows <- rows_in_period[is.finite(nee_vec[rows_in_period])]
  if (length(train_rows) < 10)
    return(list(reco = rep(NA_real_, length(rows_in_period)),
                gpp  = rep(NA_real_, length(rows_in_period))))

  night_rows <- train_rows[is_night[train_rows] & nee_vec[train_rows] > 0]
  R10 <- if (length(night_rows) > 0) {
    arr <- arrhenius_temp_scaling(temp_vec[night_rows])
    d   <- sum(arr^2, na.rm = TRUE)
    if (is.finite(d) && d > 0) sum(nee_vec[night_rows] * arr, na.rm = TRUE) / d
    else 1.0
  } else 1.0

  reco_period <- R10 * arrhenius_temp_scaling(temp_vec[rows_in_period])

  day_rows <- train_rows[is_day[train_rows] & nee_vec[train_rows] < 0]
  if (length(day_rows) > 10) {
    reco_day <- R10 * arrhenius_temp_scaling(temp_vec[day_rows])
    gpp_obs  <- reco_day - nee_vec[day_rows]
    opt <- tryCatch(
      optim(par = gpp_init, fn = gpp_ssr,
            ppfd_vec = ppfd_vec[day_rows], gpp_obs = gpp_obs, method = "BFGS"),
      error = function(e) list(par = gpp_init))
    gpp_params <- if (!is.null(opt$par)) opt$par else gpp_init
  } else {
    gpp_params <- gpp_init
  }

  gpp_raw <- light_response_gpp(gpp_params, ppfd_vec[rows_in_period])
  gpp_raw[!is_day[rows_in_period]]         <- 0
  gpp_raw[is.na(ppfd_vec[rows_in_period])] <- NA_real_

  list(reco = reco_period, gpp = gpp_raw)
}

# Internal helper: apply .partition_one_period across all regrowth periods
# and handle the fallback for period 0 (before first grazing event)
.partition_full_series <- function(nee_vec, ppfd_vec, temp_vec, regrowth,
                                    rg_ids, ppfd_night_thr, gpp_init) {
  n        <- length(nee_vec)
  is_night <- is.finite(ppfd_vec) & ppfd_vec < ppfd_night_thr
  is_day   <- is.finite(ppfd_vec) & ppfd_vec > ppfd_night_thr
  reco_out <- rep(NA_real_, n)
  gpp_out  <- rep(NA_real_, n)

  for (rg_id in rg_ids) {
    rip <- which(regrowth == rg_id)
    res <- .partition_one_period(rip, nee_vec, ppfd_vec, temp_vec,
                                  is_night, is_day, gpp_init, ppfd_night_thr)
    reco_out[rip] <- res$reco
    gpp_out[rip]  <- res$gpp
  }

  # Fallback: global fit for rows not covered by any named period
  fb <- which(is.na(reco_out) & is.finite(temp_vec))
  if (length(fb) > 0) {
    all_rows <- seq_len(n)
    res_fb <- .partition_one_period(all_rows, nee_vec, ppfd_vec, temp_vec,
                                     is_night, is_day, gpp_init, ppfd_night_thr)
    reco_out[fb] <- res_fb$reco[fb]
    gpp_out[fb]  <- res_fb$gpp[fb]
  }

  list(reco = reco_out, gpp = gpp_out)
}

# ---- Ground truth: Reco_orig and GPP_orig from NEE_orig --------------------
compute_reco_gpp_from_nee_orig <- function(flux_data, ppfd_night_thr = 10,
                                            gpp_init = c(0.08, 15, 0.5)) {
  if (!all(c("PPFD","Temp","regrowth_id","NEE_orig") %in% names(flux_data))) {
    log_msg("compute_reco_gpp_from_nee_orig: required columns missing — skipping.")
    return(flux_data)
  }
  res <- .partition_full_series(
    nee_vec  = as.numeric(flux_data$NEE_orig),
    ppfd_vec = as.numeric(flux_data$PPFD),
    temp_vec = as.numeric(flux_data$Temp),
    regrowth = flux_data$regrowth_id,
    rg_ids   = sort(unique(stats::na.omit(flux_data$regrowth_id))),
    ppfd_night_thr = ppfd_night_thr,
    gpp_init       = gpp_init
  )
  flux_data$Reco_orig <- res$reco
  flux_data$GPP_orig  <- res$gpp
  log_msg("Reco_orig and GPP_orig computed from NEE_orig (ground truth, per regrowth period).")
  flux_data
}

# ---- Per-gap: reconstruct NEE_Si_with_predicted, partition, write Si rows --
derive_reco_gpp_from_filled_nee <- function(flux_data,
                                              gap_size_cat,
                                              model_key,
                                              ppfd_night_thr = 10,
                                              gpp_init       = c(0.08, 15, 0.5)) {
  nee_pred_col  <- paste0("NEE_",  gap_size_cat, "_", model_key, "_predicted")
  reco_pred_col <- paste0("Reco_", gap_size_cat, "_", model_key, "_predicted")
  gpp_pred_col  <- paste0("GPP_",  gap_size_cat, "_", model_key, "_predicted")

  if (!nee_pred_col %in% names(flux_data)) {
    log_msg("derive_reco_gpp [", gap_size_cat, " / ", model_key,
            "]: NEE prediction column missing — skipping.")
    return(flux_data)
  }

  gap_labels <- names(flux_data)[grepl(paste0("^", gap_size_cat, "\\d+$"), names(flux_data))]
  gap_labels <- gap_labels[order(as.integer(sub(gap_size_cat, "", gap_labels)))]
  if (!length(gap_labels)) return(flux_data)

  ppfd_vec <- as.numeric(flux_data[["PPFD"]])
  temp_vec <- as.numeric(flux_data[["Temp"]])
  regrowth <- flux_data[["regrowth_id"]]
  rg_ids   <- sort(unique(stats::na.omit(regrowth)))
  nee_pred_cat <- as.numeric(flux_data[[nee_pred_col]])

  flux_data[[reco_pred_col]] <- NA_real_
  flux_data[[gpp_pred_col]]  <- NA_real_

  for (gap_lbl in gap_labels) {
    gap_rows <- which(flux_data[[gap_lbl]] %in% c(TRUE, 1))
    if (!length(gap_rows)) next

    # Step 1 — Reconstruct the complete NEE series for gap label Si:
    #   start from NEE_Si (observed outside gap, NA inside gap)
    #   fill the gap rows with model predictions from NEE_{SIZE}_{model}_predicted
    masked_col <- paste0("NEE_", gap_lbl)
    nee_complete <- if (masked_col %in% names(flux_data))
      as.numeric(flux_data[[masked_col]])          # NA at Si rows
    else
      flux_data$NEE_orig                           # fallback if column absent
    nee_complete[gap_rows] <- nee_pred_cat[gap_rows]   # fill Si with predictions

    # Step 2 — Partition the full reconstructed series per regrowth period
    res <- .partition_full_series(
      nee_vec  = nee_complete,
      ppfd_vec = ppfd_vec,
      temp_vec = temp_vec,
      regrowth = regrowth,
      rg_ids   = rg_ids,
      ppfd_night_thr = ppfd_night_thr,
      gpp_init       = gpp_init
    )

    # Step 3 — Write only Si rows into category-level prediction columns
    flux_data[[reco_pred_col]][gap_rows] <- res$reco[gap_rows]
    flux_data[[gpp_pred_col]][gap_rows]  <- res$gpp[gap_rows]
  }

  log_msg("Reco/GPP from reconstructed NEE_Si_with_predicted [",
          gap_size_cat, " / ", model_key, "].")
  flux_data
}

# =============================================================================
# SECTION 9 — MDS cross-validation  (NEE only)
# =============================================================================
# Two-pass design per gap label:
#   Pass 1 (baseline): MDS on NEE_orig — fills real NAs in the final output.
#   Pass 2 (per-gap):  MDS on NEE_{label} (NA at gap rows) — CV predictions.
# Predictions written ONLY to gap rows (not to observed rows).

run_mds_for_gap_size <- function(df_cv, gap_size_cat) {
  log_msg("MDS [", gap_size_cat, "]: starting.")
  gap_labels <- names(df_cv)[grepl(paste0("^",gap_size_cat,"\\d+$"),names(df_cv))]
  gap_labels <- gap_labels[order(as.integer(sub(gap_size_cat,"",gap_labels)))]
  if (!length(gap_labels)) {
    log_msg("MDS [",gap_size_cat,"]: no gap columns — skipping."); return(df_cv)
  }
  log_msg("MDS [",gap_size_cat,"]: ",length(gap_labels)," gap labels.")

  pred_col           <- paste0("NEE_", gap_size_cat, "_mds_predicted")
  df_cv[[pred_col]]   <- NA_real_

  # Baseline run on NEE_orig — used to fill real-NA rows
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

  # Per-gap runs
  for (gap_lbl in gap_labels) {
    masked_nee <- paste0("NEE_", gap_lbl)
    if (!masked_nee %in% names(df_cv)) next
    gap_rows <- which(df_cv[[gap_lbl]] %in% c(TRUE,1)); if (!length(gap_rows)) next

    gap_preds <- tryCatch({
      eddy_in <- prepare_reddyproc_input(df_cv, target_col = masked_nee)
      res     <- run_one_mds(eddy_in)
      res$mds_prediction[match(df_cv$timestamp, res$DateTime)]
    }, error = function(e) {
      log_msg("MDS failed for ", gap_lbl, ": ", conditionMessage(e))
      rep(NA_real_, nrow(df_cv))
    })

    df_cv[[pred_col]][gap_rows] <- gap_preds[gap_rows]
    log_msg("MDS [",gap_size_cat,"] ",gap_lbl,": filled ",length(gap_rows)," gap rows.")
  }

  log_msg("MDS [",gap_size_cat,"] complete — MAE = ",
          signif(calc_mae(df_cv$NEE_orig, df_cv[[pred_col]]),3),
          "  RMSE = ",signif(calc_rmse(df_cv$NEE_orig, df_cv[[pred_col]]),3),
          "  R² = ",  signif(calc_r2(  df_cv$NEE_orig, df_cv[[pred_col]]),3))
  df_cv
}


# =============================================================================
# SECTION 10 — Construct artificial gaps
# =============================================================================
gap_cfg <- list(
  VL = list(base=30*48, tol=250, prefix="VL"),
  L  = list(base=14*48, tol=120, prefix="L"),
  M  = list(base= 7*48, tol= 60, prefix="M"),
  S  = list(base= 3*48, tol= 10, prefix="S")
)
log_msg("Constructing artificial gaps on finite-NEE subset (VL → L → M → S) ...")
cv_VL <- create_cv_gaps_on_finite_subset(flux_data_finite, df,           gap_cfg$VL$base, gap_cfg$VL$tol, "VL")
cv_L  <- create_cv_gaps_on_finite_subset(flux_data_finite, cv_VL$data,   gap_cfg$L$base,  gap_cfg$L$tol,  "L")
cv_M  <- create_cv_gaps_on_finite_subset(flux_data_finite, cv_L$data,    gap_cfg$M$base,  gap_cfg$M$tol,  "M")
cv_S  <- create_cv_gaps_on_finite_subset(flux_data_finite, cv_M$data,    gap_cfg$S$base,  gap_cfg$S$tol,  "S")
df_cv <- cv_S$data
log_msg("Gaps ready — VL:",length(cv_VL$sizes)," L:",length(cv_L$sizes),
        " M:",length(cv_M$sizes)," S:",length(cv_S$sizes))


# =============================================================================
# SECTION 11 — Assign regrowth periods (for post-hoc Reco/GPP derivation)
# =============================================================================
if (RECO_GPP_ENABLED) {
  if ("Grazing_days_since" %in% names(df_cv)) {
    df_cv <- df_cv %>%
      dplyr::mutate(regrowth_id = build_regrowth_periods(Grazing_days_since, threshold=1.0))
    log_msg("Regrowth periods: ", dplyr::n_distinct(df_cv$regrowth_id), " unique periods.")
  } else {
    df_cv$regrowth_id <- 1L
    log_msg("Grazing_days_since absent — one global regrowth period.")
  }
  # Compute ground-truth Reco_orig / GPP_orig from NEE_orig per regrowth period
  df_cv <- compute_reco_gpp_from_nee_orig(df_cv)
}


# =============================================================================
# SECTION 12 — Run MDS cross-validation  (NEE only)
# =============================================================================
log_msg("=== Running MDS cross-validation — NEE ===")
for (cat in c("VL","L","M","S"))
  df_cv <- run_mds_for_gap_size(df_cv, cat)
log_msg("All MDS NEE gap sizes complete.")


# =============================================================================
# SECTION 13 — Derive Reco and GPP from gap-filled NEE
# =============================================================================
if (RECO_GPP_ENABLED) {
  log_msg("=== Deriving Reco and GPP from MDS gap-filled NEE ===")
  for (cat in c("VL","L","M","S"))
    df_cv <- derive_reco_gpp_from_filled_nee(df_cv, cat, model_key="mds")
  log_msg("Reco/GPP derivation complete.")
}


# =============================================================================
# SECTION 14 — Save output
# =============================================================================
dir.create(RESULTS_DIR, recursive=TRUE, showWarnings=FALSE)
saveRDS(df_cv, file.path(RESULTS_DIR, "df_cv_all_predictions.rds"))
log_msg("Saved: df_cv_all_predictions.rds")
writeLines(c(
  paste("Run finished    :", as.character(Sys.time())),
  paste("R version       :", R.version.string),
  paste("Source RDS      :", rds_name),
  paste("RESULTS_DIR     :", RESULTS_DIR),
  paste("n_rows_full     :", nrow(df_cv)),
  paste("n_rows_finite   :", nrow(flux_data_finite)),
  paste("Targets         :", if(RECO_GPP_ENABLED) "NEE (MDS), Reco+GPP (post-hoc)" else "NEE only")
), file.path(RESULTS_DIR, "run_info.txt"))
log_msg("Saved: run_info.txt — MDS script finished.")
# ================================ end MDS_CV.R ================================
