# =============================================================================
# RF_CV.R — Random Forest Cross-Validation Gap-Filling
#            NEE gap-fill with post-hoc Reco and GPP derivation
# =============================================================================
#
# PURPOSE
#   Trains a Random Forest model (R package {ranger}) separately for each
#   artificial gap in each of four gap-size categories (S, M, L, VL).
#   For each gap, the model is trained on all observations outside that
#   gap where NEE_orig is finite, then used to predict NEE at the gap rows.
#   Reco and GPP are subsequently derived from the gap-filled NEE series
#   through flux partitioning (see Section 7).
#
# WHY NEE IS GAP-FILLED FIRST
#   Reco and GPP are diagnostic quantities derived from NEE via a process
#   model; they are not independently measured.  Gap-filling NEE first and
#   then partitioning preserves the carbon-balance identity NEE = Reco - GPP
#   at every gap row by construction.  Estimating all three variables
#   independently risks violating this constraint unless explicit mass-
#   balance corrections are applied.
#
# INPUTS (passed from the calling run script)
#   df          — prepared data frame (data/data_prepared/{SITE}.rds)
#   predictors  — character vector of predictor column names
#   RESULTS_DIR — output directory path
#   rds_name    — source file name (written to run_info.txt)
#
# OUTPUTS  (written to RESULTS_DIR)
#   df_cv_all_predictions.rds
#     The input data frame with 12 new columns:
#       NEE_{S|M|L|VL}_rf_predicted
#       Reco_{S|M|L|VL}_rf_predicted
#       GPP_{S|M|L|VL}_rf_predicted
#
#   variable_importance/
#     rf_vi_{SIZE}_{LABEL}.rds           — per-gap importance list
#     rf_variable_importance_{SIZE}_all.rds / .csv — stacked table
#
#   progress.log | progress_log.rds | run_info.txt
#
# SECTIONS
#    1  Reproducibility           — random seed and thread control
#    2  Required packages         — package availability check
#    3  Progress logger           — timestamped log to console and file
#    4  Helper functions          — gap partitioner, gap creator, metrics
#    5  Phytomass Index           — per-gap PI computation (managed sites)
#    6  Regrowth period detection — grazing-event based period assignment
#    7  Flux partitioning helpers — Reco/GPP models and partitioning logic
#    8  RF cross-validation       — run_rf_for_gap_size() definition
#    9  Data preparation          — filter finite NEE; assign regrowth IDs
#   10  Construct artificial gaps — build S/M/L/VL gap matrices
#   11  Phytomass Index           — add PI columns if enabled
#   12  Run RF cross-validation   — loop over gap sizes
#   13  Derive Reco and GPP       — post-hoc partitioning from gap-filled NEE
#   14  Save outputs              — write RDS and run_info.txt
#
# NOTE ON FEATURE SCALING
#   Random Forest splits on thresholds are invariant to monotone
#   transformations of predictors.  Min-max scaling is therefore not
#   required and is deliberately omitted.  See Hastie et al. (2009) §15.2.
#
# =============================================================================


# =============================================================================
# SECTION 1 — Reproducibility
# =============================================================================
# A fixed random seed ensures that identical results are obtained across
# runs on the same machine.  Thread counts for OpenMP, MKL, and OpenBLAS
# are fixed at 1 to prevent non-determinism from parallel floating-point
# reductions, which can produce slightly different results depending on
# thread scheduling.

set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")


# =============================================================================
# SECTION 2 — Required packages
# =============================================================================
# All required packages are checked for availability before any computation
# begins.  A single informative error is raised listing all missing packages,
# so the user can install them in one step rather than discovering them
# one at a time.

needed <- c("dplyr", "tibble", "tidyr", "ggplot2", "glue",
            "lubridate", "rlang", "purrr", "ranger")
miss <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))


# =============================================================================
# SECTION 3 — Progress logger
# =============================================================================
# The log_msg() function writes timestamped messages to three locations
# simultaneously: (1) the R console for interactive monitoring, (2) a plain-
# text progress.log file for post-hoc inspection, and (3) an RDS list for
# programmatic access.  This redundancy ensures progress is recoverable
# even if the console output is lost in a batch job.
# RESULTS_DIR falls back to a temporary directory if not passed by the
# calling script, making the file sourceable in isolation for testing.

if (!exists("RESULTS_DIR", inherits = TRUE) || is.null(RESULTS_DIR))
  RESULTS_DIR <- file.path(tempdir(), "rf_cv_fallback")
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
log_msg("RF (NEE gap-fill + post-hoc Reco/GPP) started.  RESULTS_DIR = ", RESULTS_DIR)


# =============================================================================
# SECTION 4 — Helper functions
# =============================================================================
# Four utilities shared by the gap construction and evaluation steps.
# They are defined once here and used throughout the script.

# 4a  Block partitioner — partition_timeseries_into_blocks()
# ---------------------------------------------------------------------------
# Divides a time series of n rows into k contiguous blocks each of
# approximately base_size rows.  Remainder rows are distributed across
# the last min(k, ceil(remainder/tol)) blocks, subject to the constraint
# that no block deviates from base_size by more than tol rows.  This
# ensures that all artificial gaps are of comparable duration while
# exactly covering the entire dataset.  The function returns both the
# block sizes and a list of integer index vectors (one per block).
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

# 4b  Artificial gap creator — create_cv_gaps()
# ---------------------------------------------------------------------------
# Wraps partition_timeseries_into_blocks() and builds two matrices:
#   fm — Boolean flag matrix (n × k): fm[i,j] = TRUE if row i is inside gap j
#   nm — Masked NEE matrix (n × k): NEE_orig with gap rows set to NA
# Both are appended to the data frame as columns named {prefix}N and
# NEE_{prefix}N respectively.  The flag columns identify which rows belong
# to each gap; the masked columns provide the training target for that gap
# (observed NEE outside the gap, NA inside it).
create_cv_gaps <- function(flux_data, base_size, tol, prefix) {
  n <- nrow(flux_data)
  bl <- partition_timeseries_into_blocks(n, as.integer(round(base_size)),
                                          as.integer(round(tol)))
  fn <- paste0(prefix, seq_along(bl$indices)); nn <- paste0("NEE_", fn)
  fm <- matrix(FALSE,    nrow = n, ncol = length(fn), dimnames = list(NULL, fn))
  nm <- matrix(NA_real_, nrow = n, ncol = length(fn), dimnames = list(NULL, nn))
  for (i in seq_along(bl$indices)) {
    fm[bl$indices[[i]], i] <- TRUE
    nm[, i] <- flux_data$NEE_orig; nm[bl$indices[[i]], i] <- NA_real_
  }
  list(data  = dplyr::bind_cols(flux_data,
                                 tibble::as_tibble(as.data.frame(fm)),
                                 tibble::as_tibble(as.data.frame(nm))),
       sizes = bl$sizes)
}

# 4c  Metric functions — calc_mae(), calc_rmse(), calc_r2()
# ---------------------------------------------------------------------------
# All three functions first restrict to row pairs where both observed (o)
# and predicted (p) values are finite (.vp() guard).  This is essential
# because gaps in Reco_orig and GPP_orig can produce NA ground-truth values
# outside the originally intended gap rows.  Returning NA rather than a
# biased metric when no valid pairs exist prevents silent propagation of
# NaN values into summary tables.
.vp <- function(o, p) is.finite(as.numeric(o)) & is.finite(as.numeric(p))
calc_mae  <- function(o, p) { ok<-.vp(o,p); if(!any(ok)) return(NA_real_); mean(abs(as.numeric(p)[ok]-as.numeric(o)[ok])) }
calc_rmse <- function(o, p) { ok<-.vp(o,p); if(!any(ok)) return(NA_real_); sqrt(mean((as.numeric(p)[ok]-as.numeric(o)[ok])^2)) }
calc_r2   <- function(o, p) { ok<-.vp(o,p); if(sum(ok)<2) return(NA_real_); ov<-as.numeric(o)[ok]; pv<-as.numeric(p)[ok]; den<-sum((ov-mean(ov))^2)*sum((pv-mean(pv))^2); if(den<=0) return(NA_real_); (sum((ov-mean(ov))*(pv-mean(pv)))^2)/den }

# 4d  Season label helper — timestamp_to_season()
# ---------------------------------------------------------------------------
# Assigns a meteorological season label to each timestamp.  The convention
# used here (Nov-Jan = Winter, Feb-Apr = Spring, May-Jul = Summer,
# Aug-Oct = Autumn) is based on the Irish meteorological calendar.
timestamp_to_season <- function(x) {
  m <- as.integer(format(x, "%m"))
  dplyr::case_when(m %in% c(11,12,1) ~ "Winter", m %in% c(2,3,4) ~ "Spring",
                   m %in% c(5,6,7)   ~ "Summer", m %in% c(8,9,10) ~ "Autumn",
                   TRUE ~ NA_character_)
}


# =============================================================================
# SECTION 5 — Phytomass Index  (managed grassland only)
# =============================================================================
# The Phytomass Index (PI) is a derived predictor that estimates the state
# of sward regrowth at the time of each artificial gap.  It is computed
# from night-time and daytime NEE observations OUTSIDE each gap over a
# rolling window of window_days (default 21 days), capturing the net carbon
# balance of the canopy in the period surrounding the gap.
#
# PI is activated when "Grazing_days_since" is present in the predictors
# vector, because regrowth periods — which define the temporal structure
# of PI — are delineated by grazing events.  One PI_{gap_label} column is
# added per gap label, so the model sees a gap-specific regrowth state
# rather than a single global value.
#
# The add_phytomass_index() function implements the following steps:
#   1. Aggregate daily night NEE (PPFD < night_ppfd) and day NEE
#      (PPFD > day_ppfd) across the full time series (dt) and within
#      gap rows only (gd) for each gap label.
#   2. Compute the outside-gap daily averages as dt - gd (subtracting the
#      gap contribution from the full totals preserves independence).
#   3. Apply a centred rolling sum (window_days) to both day and night
#      outside-gap averages.
#   4. PI_raw = rolling_mean_night_NEE - rolling_mean_day_NEE.
#      Night NEE is positive (respiration) and day NEE is negative
#      (photosynthesis), so PI_raw increases as the sward grows.
#   5. Normalise PI to [0, 1] by dividing by the maximum observed PI_raw.
#   6. Set PI to NA at gap rows (no information available inside the gap)
#      and linearly interpolate through those NA values.
# PI_ENABLED is derived from the predictor list; it can also be pre-set
# by the calling run script to override the default behaviour.

PI_ENABLED <- "Grazing_days_since" %in% get0("predictors", ifnotfound = character(0))
if (!PI_ENABLED) log_msg("Phytomass Index disabled (Grazing_days_since not in predictors).")

.find_ppfd_col <- function(df) {
  p <- c("PPFD","PPFD_IN","PPFD_IN_F","PPFD_F","PPFD_1_1_1_gapfilled","PPFD_1_1_1")
  h <- p[p %in% names(df)]; if (length(h)) return(h[1])
  h2 <- grep("PPFD", names(df), ignore.case = TRUE, value = TRUE)
  if (length(h2)) return(h2[1]); NULL
}

add_phytomass_index <- function(flux_data, gap_prefixes = c("VL","L","M","S"),
                                 ppfd_col = NULL, night_ppfd = 1, day_ppfd = 700,
                                 window_days = 21L) {
  if (!requireNamespace("zoo", quietly = TRUE))
    stop("Package 'zoo' required for PI.  Run: install.packages('zoo')")
  if (!inherits(flux_data$timestamp, "POSIXt"))
    flux_data$timestamp <- suppressWarnings(
      lubridate::parse_date_time(flux_data$timestamp,
        orders = c("Ymd HMS","Ymd HM","Ymd","Y/m/d HMS","Y/m/d")))
  if (is.null(ppfd_col)) ppfd_col <- .find_ppfd_col(flux_data)
  if (is.null(ppfd_col) || !ppfd_col %in% names(flux_data)) {
    log_msg("PI: no PPFD column — skipping."); return(flux_data)
  }
  window_days <- as.integer(window_days)
  if (!is.finite(window_days) || window_days < 1) window_days <- 21L
  if (window_days %% 2 == 0) window_days <- window_days + 1L

  flux_data$.cd <- as.Date(flux_data$timestamp)
  nv <- as.numeric(flux_data$NEE_orig)
  pv <- suppressWarnings(as.numeric(flux_data[[ppfd_col]]))
  ok <- is.finite(nv)&is.finite(pv); isn <- ok&pv<night_ppfd; isd <- ok&pv>day_ppfd

  dt <- tibble::tibble(date=flux_data$.cd,
                       nn=dplyr::if_else(isn,nv,0), cn=as.integer(isn),
                       nd=dplyr::if_else(isd,nv,0), cd=as.integer(isd)) %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(nn=sum(nn,na.rm=T), cn=sum(cn,na.rm=T),
                     nd=sum(nd,na.rm=T), cd=sum(cd,na.rm=T), .groups="drop") %>%
    dplyr::arrange(date)
  rol <- function(x) zoo::rollapply(x, window_days, function(v) sum(v,na.rm=T),
                                    align="center", fill=NA_real_)

  build_pi <- function(gl) {
    gi <- which(flux_data[[gl]] %in% c(TRUE,1)); if (!length(gi)) return(NULL)
    gd <- tibble::tibble(date=flux_data$.cd[gi],
                         nn=dplyr::if_else(isn[gi],nv[gi],0), cn=as.integer(isn[gi]),
                         nd=dplyr::if_else(isd[gi],nv[gi],0), cd=as.integer(isd[gi])) %>%
      dplyr::group_by(date) %>%
      dplyr::summarise(nn=sum(nn,na.rm=T), cn=sum(cn,na.rm=T),
                       nd=sum(nd,na.rm=T), cd=sum(cd,na.rm=T), .groups="drop")
    dd <- dplyr::left_join(dt, gd, by="date", suffix=c("","_g")) %>%
      dplyr::mutate(across(ends_with("_g"), ~dplyr::coalesce(.,0)),
                    an=pmax(nn-nn_g,0), acn=pmax(cn-cn_g,0L),
                    ad=pmax(nd-nd_g,0), acd=pmax(cd-cd_g,0L)) %>%
      dplyr::arrange(date)
    mn <- dplyr::if_else(rol(dd$acn)>0, rol(dd$an)/rol(dd$acn), NA_real_)
    md <- dplyr::if_else(rol(dd$acd)>0, rol(dd$ad)/rol(dd$acd), NA_real_)
    PR <- dplyr::if_else(is.finite(mn)&is.finite(md), mn-md, NA_real_)
    mx <- suppressWarnings(max(PR, na.rm=T))
    PN <- if (!is.finite(mx)||mx<=0) rep(NA_real_,length(PR)) else pmax(0,pmin(PR/mx,1))
    ph <- PN[match(flux_data$.cd, dd$date)]; ph[gi] <- NA_real_
    ord <- order(flux_data$timestamp); po <- ph[ord]; xo <- as.numeric(flux_data$timestamp[ord])
    rl  <- c(if(is.na(po[1]))2 else 1, if(is.na(po[length(po)]))2 else 1)
    if (sum(is.finite(po))>=2) po <- zoo::na.approx(po, x=xo, na.rm=FALSE, rule=rl)
    ph[ord] <- po; pf <- PN[match(flux_data$.cd, dd$date)]; pf[gi] <- ph[gi]; pf
  }

  for (pfx in gap_prefixes) {
    gls <- names(flux_data)[grepl(paste0("^",pfx,"\\d+$"), names(flux_data))]
    gls <- gls[order(as.integer(sub(pfx,"",gls)))]; if (!length(gls)) next
    log_msg("PI: computing ", length(gls), " columns for '", pfx, "'.")
    for (gl in gls) { pv2 <- build_pi(gl); if (!is.null(pv2)) flux_data[[paste0("PI_",gl)]] <- pv2 }
  }
  flux_data$.cd <- NULL; flux_data
}


# =============================================================================
# SECTION 6 — Regrowth period detection
# =============================================================================
# Assigns an integer period ID to every row, where each new period begins
# at the first observation after a grazing event (Grazing_days_since == 0
# following a row where days_since was > threshold or missing).
# The regrowth_id is used by the flux partitioning functions in Section 7
# to fit separate Reco and GPP parameters for each regrowth cycle, because
# post-grazing canopy structure (LAI, biomass) strongly influences both
# light-use efficiency and the temperature sensitivity of respiration.

build_regrowth_periods <- function(days_since, threshold = 1.0) {
  x       <- suppressWarnings(as.numeric(days_since))
  is_zero <- is.finite(x) & (x == 0); x_prev <- dplyr::lag(x)
  cumsum(ifelse(is_zero & (is.na(x_prev)|!is.finite(x_prev)|(x_prev>threshold)), 1L, 0L))
}


# =============================================================================
# SECTION 7 — Flux partitioning helpers
# =============================================================================
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
# SECTION 8 — Random Forest cross-validation  (NEE only)
# =============================================================================
# run_rf_for_gap_size() iterates over every gap label within one gap-size
# category (e.g. all "M1", "M2", … gaps).  For each label:
#
#   1. Identify training rows: rows where NEE_{gap_label} is finite
#      (i.e. all rows OUTSIDE the current gap).
#   2. If PI is enabled, append PI_{gap_label} to the feature set so that
#      the model has a gap-specific regrowth index as a predictor.
#   3. Train ranger() on the training rows with impurity importance.
#   4. Predict NEE at the gap rows and write into NEE_{SIZE}_rf_predicted.
#   5. Save per-gap variable importance as rf_vi_{SIZE}_{LABEL}.rds.
#
# After all gap labels are processed, the per-gap importance tables are
# stacked with bind_rows() and saved as the aggregate file
# rf_variable_importance_{SIZE}_all.rds / .csv.
#
# Arguments:
#   flux_data    — data frame with gap flag and masked NEE columns
#   gap_size_cat — "S" | "M" | "L" | "VL"
#   feature_cols — predictor column names
#   random_seed  — passed to ranger(); default 42 for reproducibility
#   importance_type — "impurity" (Gini) or "none" to skip VI saving
# Returns:
#   flux_data with column  NEE_{gap_size_cat}_rf_predicted  added

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

  prediction_col             <- paste0("NEE_", gap_size_cat, "_rf_predicted")
  flux_data[[prediction_col]] <- NA_real_

  vi_dir <- file.path(RESULTS_DIR, "variable_importance")
  dir.create(vi_dir, recursive = TRUE, showWarnings = FALSE)
  vi_list <- list()

  for (gap_lbl in gap_labels) {
    masked_nee <- paste0("NEE_", gap_lbl)
    if (!masked_nee %in% names(flux_data)) next

    feats <- feature_cols
    if (PI_ENABLED) {
      pi_col <- paste0("PI_", gap_lbl)
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
                                           gap_size = gap_size_cat,
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
# SECTION 9 — Data preparation
# =============================================================================
# The model is trained only on rows where NEE_orig is finite.  Rows with
# missing NEE (genuine measurement gaps) are excluded from training but
# are retained in the output data frame so that predictions can be written
# back to the correct positions.  RECO_GPP_ENABLED is checked once here;
# flux partitioning is skipped entirely if PPFD or Temp is absent.

flux_data <- df %>% dplyr::filter(is.finite(NEE_orig))
log_msg("Rows with finite NEE_orig retained: ", nrow(flux_data))

RECO_GPP_ENABLED <- all(c("PPFD", "Temp") %in% names(flux_data))
if (!RECO_GPP_ENABLED)
  log_msg("WARNING: PPFD or Temp missing — Reco/GPP derivation will be skipped.")

if (RECO_GPP_ENABLED) {
  if ("Grazing_days_since" %in% names(flux_data)) {
    flux_data <- flux_data %>%
      dplyr::mutate(regrowth_id = build_regrowth_periods(Grazing_days_since, threshold = 1.0))
    log_msg("Regrowth periods: ", dplyr::n_distinct(flux_data$regrowth_id), " unique periods.")
  } else {
    flux_data$regrowth_id <- 1L
    log_msg("Grazing_days_since absent — using one global regrowth period.")
  }
  # Compute ground-truth Reco_orig / GPP_orig from NEE_orig per regrowth period
  flux_data <- compute_reco_gpp_from_nee_orig(flux_data)
}


# =============================================================================
# SECTION 10 — Construct artificial gaps
# =============================================================================
# Four gap-size categories are constructed sequentially on the same data
# frame, so that all gap matrices coexist in a single object:
#   VL  — very long (~30 days = 30 × 48 half-hours), tolerance ±250 rows
#   L   — long     (~14 days = 14 × 48 half-hours), tolerance ±120 rows
#   M   — medium   (~7  days =  7 × 48 half-hours), tolerance ±60  rows
#   S   — short    (~3  days =  3 × 48 half-hours), tolerance ±10  rows
# The sequential construction means that M gaps are created from the data
# already containing VL and L gap columns — they all share the same row
# indices, ensuring gap labels are temporally consistent across categories.

gap_cfg <- list(
  VL = list(base = 30*48, tol = 250, prefix = "VL"),
  L  = list(base = 14*48, tol = 120, prefix = "L"),
  M  = list(base =  7*48, tol =  60, prefix = "M"),
  S  = list(base =  3*48, tol =  10, prefix = "S")
)
log_msg("Constructing artificial gaps (VL → L → M → S) ...")
cv_VL <- create_cv_gaps(flux_data,   gap_cfg$VL$base, gap_cfg$VL$tol, "VL")
cv_L  <- create_cv_gaps(cv_VL$data,  gap_cfg$L$base,  gap_cfg$L$tol,  "L")
cv_M  <- create_cv_gaps(cv_L$data,   gap_cfg$M$base,  gap_cfg$M$tol,  "M")
cv_S  <- create_cv_gaps(cv_M$data,   gap_cfg$S$base,  gap_cfg$S$tol,  "S")
flux_data <- cv_S$data
log_msg("Gaps ready — VL:", length(cv_VL$sizes), " L:", length(cv_L$sizes),
        " M:", length(cv_M$sizes), " S:", length(cv_S$sizes))


# =============================================================================
# SECTION 11 — Phytomass Index
# =============================================================================
# add_phytomass_index() is called AFTER all gap matrices are constructed
# so that it can compute a separate PI value for each gap label.  The
# window_days parameter (21 days = 21-day centred rolling window) controls
# the temporal smoothing of the night/day NEE averages used to estimate
# regrowth state.  Odd values are enforced for symmetric centring.

if (PI_ENABLED) {
  flux_data <- add_phytomass_index(flux_data, gap_prefixes = c("VL","L","M","S"),
                                   window_days = 21L)
  log_msg("PI columns added.")
}


# =============================================================================
# SECTION 12 — Run RF cross-validation  (NEE only)
# =============================================================================
# The RF is trained and evaluated for NEE across all four gap-size categories.
# Gap sizes are processed largest-first (VL → L → M → S) so that longer-
# gap training sets are available early; order does not affect correctness
# because each gap label is an independent training/prediction split.

log_msg("=== Running RF cross-validation — NEE ===")
for (cat in c("VL","L","M","S"))
  flux_data <- run_rf_for_gap_size(flux_data, cat, predictors)
log_msg("All RF NEE gap sizes complete.")


# =============================================================================
# SECTION 13 — Derive Reco and GPP from gap-filled NEE
# =============================================================================
# For each gap size category, derive_reco_gpp_from_filled_nee() reconstructs
# a complete NEE time series for each gap label by combining observed NEE
# outside the gap with RF predictions inside it.  The full reconstructed
# series is then partitioned per regrowth period to produce Reco and GPP
# at the gap rows.  Only the gap rows are written to the output columns;
# the observed-period partitioning is not stored.

if (RECO_GPP_ENABLED) {
  log_msg("=== Deriving Reco and GPP from RF gap-filled NEE ===")
  for (cat in c("VL","L","M","S"))
    flux_data <- derive_reco_gpp_from_filled_nee(flux_data, cat, model_key = "rf")
  log_msg("Reco/GPP derivation complete.")
}


# =============================================================================
# SECTION 14 — Save outputs
# =============================================================================
# The complete data frame (original columns + all 12 prediction columns +
# gap flag and masked NEE columns) is saved as df_cv_all_predictions.rds.
# The run_info.txt file records key metadata for reproducibility and
# provenance tracking.

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(flux_data, file.path(RESULTS_DIR, "df_cv_all_predictions.rds"))
log_msg("Saved: df_cv_all_predictions.rds")

writeLines(c(
  paste("Run finished    :", as.character(Sys.time())),
  paste("R version       :", R.version.string),
  paste("Platform        :", R.version$platform),
  paste("Source RDS      :", rds_name),
  paste("RESULTS_DIR     :", RESULTS_DIR),
  paste("n_rows          :", nrow(flux_data)),
  paste("Targets         :", if (RECO_GPP_ENABLED) "NEE (RF), Reco+GPP (post-hoc)" else "NEE only"),
  paste("Predictors      :", paste(predictors, collapse = ", "))
), file.path(RESULTS_DIR, "run_info.txt"))
log_msg("Saved: run_info.txt — RF script finished.")

# ============================== end RF_CV.R ===================================
