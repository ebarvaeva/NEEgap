# =============================================================================
# miniRECgap_CV.R — Process-Based Gap-Filling: Reco and GPP from First Principles
# =============================================================================
#
# PURPOSE
#   Fills artificial gaps in NEE using a process-based model that derives
#   ecosystem respiration (Reco) and gross primary production (GPP) from
#   fitted semi-empirical relationships, then computes:
#
#     NEE_pred = Reco_model - GPP_model
#
#   Parameters are estimated separately for each artificial gap label and
#   each vegetation regrowth period using only the observations OUTSIDE
#   the gap.  This preserves the cross-validation integrity of the
#   evaluation: no information from inside the gap influences the fitted
#   parameters used to fill it.
#
# PROCESS MODELS
#   Reco — Lloyd-Taylor temperature response (Arrhenius form):
#     Reco(T) = R10 × exp[309 × (1/(283.2-230) - 1/(T_K-230))]
#     R10 is estimated by ordinary least squares from nighttime training
#     rows (PPFD < 10 µmol m⁻² s⁻¹) where Reco ≈ NEE > 0.
#
#   GPP  — Thornley non-rectangular hyperbola:
#     GPP(Q) = [αQ + GPPmax - sqrt((αQ+GPPmax)² - 4θαQGPPmax)] / (2θ)
#     Parameters (α, GPPmax, θ) are fitted by BFGS optimisation from
#     daytime training rows (PPFD > 10 µmol m⁻² s⁻¹) where GPP ≈ Reco - NEE.
#     GPP is set to 0 at night.
#
# REGROWTH PERIODS
#   If Grazing_days_since is present in the data, regrowth periods are
#   delineated by grazing events and parameters are fitted per period.
#   A single global period is used otherwise.  This is critical for
#   managed grasslands where post-grazing canopy structure strongly
#   influences both α and R10.
#
# HOW RECO AND GPP ARE STORED
#   Unlike RF/MLP/XGBoost where Reco and GPP are derived post-hoc from
#   the gap-filled NEE, miniRECgap computes Reco and GPP directly as
#   intermediate model outputs.  They are stored in the prediction columns
#   at the gap rows alongside NEE_pred.  This ensures NEE = Reco - GPP
#   holds by construction at every gap row, consistent with the other models.
#
# INPUTS (passed from the calling run script)
#   df, RESULTS_DIR, rds_name
#   Note: predictors is not used — miniRECgap has no machine-learning feature
#   set; it uses only Temp, PPFD, and the regrowth period structure.
#
# OUTPUTS  (written to RESULTS_DIR)
#   df_cv_all_predictions.rds
#     Original data frame augmented with:
#       NEE_{S|M|L|VL}_minirec_predicted
#       Reco_{S|M|L|VL}_minirec_predicted
#       GPP_{S|M|L|VL}_minirec_predicted
#
#   progress.log | run_info.txt
#
# SECTIONS
#    1  Reproducibility            — random seed
#    2  Required packages          — package check
#    3  Progress logger            — timestamped log
#    4  Helper functions           — gap partitioner, gap creator, metrics
#    5  Regrowth period detection  — grazing-event based period assignment
#    6  Process-based flux models  — Arrhenius Reco, light-response GPP,
#                                    fitting procedure, ground-truth computation
#    7  miniRECgap cross-validation — run_minirec_for_gap_size() definition
#    8  Data preparation
#    9  Assign regrowth periods
#   10  Construct artificial gaps
#   11  Run miniRECgap cross-validation
#   12  Save outputs
#
# =============================================================================


# =============================================================================
# SECTION 1 — Reproducibility
# =============================================================================
# Fixed random seed for reproducibility.  Thread counts are fixed at 1
# for the same reason as in RF_CV.R (deterministic floating-point reduction).

set.seed(42)
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")


# =============================================================================
# SECTION 2 — Required packages
# =============================================================================
# All packages are checked before any computation begins.

needed <- c("dplyr","tibble","tidyr","lubridate","purrr","glue")
miss   <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))


# =============================================================================
# SECTION 3 — Progress logger
# =============================================================================
# Timestamped logging to console and file.  See RF_CV.R Section 3.

if (!exists("RESULTS_DIR", inherits = TRUE) || is.null(RESULTS_DIR))
  RESULTS_DIR <- file.path(tempdir(), "minirec_cv_fallback")
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
log_msg("miniRECgap (NEE + Reco/GPP captured as intermediates) started.  RESULTS_DIR = ", RESULTS_DIR)


# =============================================================================
# SECTION 4 — Helper functions
# =============================================================================
# Gap construction and metric computation utilities.
# See RF_CV.R Section 4 for detailed documentation.

# 4a  Block partitioner
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

# 4b  Artificial gap creator
create_cv_gaps <- function(flux_data, base_size, tol, prefix) {
  n      <- nrow(flux_data)
  blocks <- partition_timeseries_into_blocks(n, as.integer(round(base_size)),
                                              as.integer(round(tol)))
  fn <- paste0(prefix, seq_along(blocks$indices)); nn <- paste0("NEE_", fn)
  fm <- matrix(FALSE,    nrow=n, ncol=length(fn), dimnames=list(NULL,fn))
  nm <- matrix(NA_real_, nrow=n, ncol=length(fn), dimnames=list(NULL,nn))
  for (i in seq_along(blocks$indices)) {
    fm[blocks$indices[[i]], i] <- TRUE
    nm[, i] <- flux_data$NEE_orig; nm[blocks$indices[[i]], i] <- NA_real_
  }
  list(data  = dplyr::bind_cols(flux_data,
                                 tibble::as_tibble(as.data.frame(fm)),
                                 tibble::as_tibble(as.data.frame(nm))),
       sizes = blocks$sizes)
}

# 4c  Metrics
.vp <- function(o,p) is.finite(as.numeric(o))&is.finite(as.numeric(p))
calc_mae  <- function(o,p){ok<-.vp(o,p); if(!any(ok)) return(NA_real_); mean(abs(as.numeric(p)[ok]-as.numeric(o)[ok]))}
calc_rmse <- function(o,p){ok<-.vp(o,p); if(!any(ok)) return(NA_real_); sqrt(mean((as.numeric(p)[ok]-as.numeric(o)[ok])^2))}
calc_r2   <- function(o,p){ok<-.vp(o,p); if(sum(ok)<2) return(NA_real_); ov<-as.numeric(o)[ok]; pv<-as.numeric(p)[ok]; den<-sum((ov-mean(ov))^2)*sum((pv-mean(pv))^2); if(den<=0) return(NA_real_); (sum((ov-mean(ov))*(pv-mean(pv)))^2)/den}


# =============================================================================
# SECTION 5 — Regrowth period detection
# =============================================================================
# assign_regrowth_periods() labels each row with an integer regrowth period
# ID.  A new period begins at the first row where days_since == 0 following
# a row where it was above the threshold (default 1 day).  This delineates
# distinct post-grazing regrowth cycles, within which the process-model
# parameters (R10, α, GPPmax, θ) are assumed to be stationary.  Using
# separate parameters per regrowth period reflects the strong influence of
# canopy age on both light-use efficiency and temperature sensitivity
# of respiration.

assign_regrowth_periods <- function(days_since, threshold = 1.0) {
  x       <- suppressWarnings(as.numeric(days_since))
  is_zero <- is.finite(x) & (x == 0); x_prev <- dplyr::lag(x)
  cumsum(ifelse(is_zero & (is.na(x_prev)|!is.finite(x_prev)|(x_prev>threshold)), 1L, 0L))
}


# =============================================================================
# SECTION 6 — Process-based flux models
# =============================================================================
# Six functions implementing the Reco and GPP process models and their
# fitting procedure:
#
# arrhenius_temperature_response(temp_celsius)
#   Lloyd-Taylor (1994) Arrhenius function for ecosystem respiration:
#   A(T) = exp[309 × (1/(283.2-230) - 1/(T_K-230))]
#   where T_K = temp + 273.15.  Reco(T) = R10 × A(T).
#
# light_response_gpp(params, ppfd_vec)
#   Thornley non-rectangular hyperbola (Thornley & Johnson 1990):
#   GPP(Q) = [αQ + GPPmax - sqrt((αQ+GPPmax)² - 4θαQGPPmax)] / (2θ)
#   params = c(α, GPPmax, θ):
#     α      — initial slope (apparent quantum yield)
#     GPPmax — maximum GPP (light saturation plateau)
#     θ      — convexity parameter (0 = rectangular, 1 = Blackman)
#
# gpp_residual_sum_squares(params, ppfd_vec, gpp_obs)
#   Objective function for BFGS optimisation of the GPP light-response.
#
# fit_period_parameters(rows, nee, ppfd, temp, is_night, is_day, gpp_init)
#   Fits R10 by OLS from nighttime training rows (positive NEE, PPFD < 10),
#   then fits (α, GPPmax, θ) by BFGS from daytime training rows (NEE < 0).
#   Returns a list with R10 and gpp_params for application at all rows.
#
# compute_reco_gpp_from_nee_orig(flux_data, ...)
#   Applies fit_period_parameters() to the observed NEE series per regrowth
#   period.  Produces Reco_orig and GPP_orig as ground-truth references.
#
# run_minirec_for_gap_size() in Section 7 uses fit_period_parameters()
# on the OUTSIDE-GAP training rows for each gap label, then applies the
# fitted parameters to compute Reco and GPP at the gap rows.

# Lloyd-Taylor Arrhenius temperature response:
#   A(T) = exp(309 × (1/53.2 − 1/(T_air + 46.02)))
#   Reco = R10 × A(T),  R10 estimated by OLS from nighttime NEE > 0
#
# Non-rectangular hyperbolic light response (Thornley 1998):
#   GPP = (α·Q + GPPmax − √((α·Q+GPPmax)² − 4·θ·α·Q·GPPmax)) / (2θ)
#   Q = PPFD,  parameters (α, GPPmax, θ) optimised by BFGS
#   GPP = 0 at night (PPFD < 10 µmol m⁻² s⁻¹)

arrhenius_temperature_response <- function(temp_celsius) {
  exp(309 * ((1 / (283.2 - 230)) - (1 / ((temp_celsius + 273.2) - 230))))
}

light_response_gpp <- function(params, ppfd_vec) {
  a <- params[1]; b <- params[2]; c <- params[3]
  num  <- a * ppfd_vec + b
  disc <- pmax(num^2 - 4 * c * (ppfd_vec * a * b), 0)
  (num - sqrt(disc)) / (2 * c)
}

gpp_residual_sum_squares <- function(params, ppfd_vec, gpp_obs) {
  sum((light_response_gpp(params, ppfd_vec) - gpp_obs)^2)
}

# Aliases so shared helpers below use same names as the other model scripts
arrhenius_temp_scaling <- arrhenius_temperature_response
gpp_ssr <- gpp_residual_sum_squares

# Shared internal helpers: partition a complete NEE series per regrowth period
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
    if (is.finite(d) && d > 0) sum(nee_vec[night_rows] * arr, na.rm = TRUE) / d else 1.0
  } else 1.0
  reco_period <- R10 * arrhenius_temp_scaling(temp_vec[rows_in_period])
  day_rows <- train_rows[is_day[train_rows] & nee_vec[train_rows] < 0]
  if (length(day_rows) > 10) {
    reco_day <- R10 * arrhenius_temp_scaling(temp_vec[day_rows])
    gpp_obs  <- reco_day - nee_vec[day_rows]
    opt <- tryCatch(optim(par=gpp_init, fn=gpp_ssr, ppfd_vec=ppfd_vec[day_rows],
                          gpp_obs=gpp_obs, method="BFGS"),
                    error=function(e) list(par=gpp_init))
    gpp_params <- if (!is.null(opt$par)) opt$par else gpp_init
  } else gpp_params <- gpp_init
  gpp_raw <- light_response_gpp(gpp_params, ppfd_vec[rows_in_period])
  gpp_raw[!is_day[rows_in_period]]         <- 0
  gpp_raw[is.na(ppfd_vec[rows_in_period])] <- NA_real_
  list(reco = reco_period, gpp = gpp_raw)
}

.partition_full_series <- function(nee_vec, ppfd_vec, temp_vec, regrowth,
                                    rg_ids, ppfd_night_thr, gpp_init) {
  n <- length(nee_vec)
  is_night <- is.finite(ppfd_vec) & ppfd_vec < ppfd_night_thr
  is_day   <- is.finite(ppfd_vec) & ppfd_vec > ppfd_night_thr
  reco_out <- rep(NA_real_, n); gpp_out <- rep(NA_real_, n)
  for (rg_id in rg_ids) {
    rip <- which(regrowth == rg_id)
    res <- .partition_one_period(rip, nee_vec, ppfd_vec, temp_vec,
                                  is_night, is_day, gpp_init, ppfd_night_thr)
    reco_out[rip] <- res$reco; gpp_out[rip] <- res$gpp
  }
  fb <- which(is.na(reco_out) & is.finite(temp_vec))
  if (length(fb) > 0) {
    res_fb <- .partition_one_period(seq_len(n), nee_vec, ppfd_vec, temp_vec,
                                     is_night, is_day, gpp_init, ppfd_night_thr)
    reco_out[fb] <- res_fb$reco[fb]; gpp_out[fb] <- res_fb$gpp[fb]
  }
  list(reco = reco_out, gpp = gpp_out)
}

# Ground truth: Reco_orig and GPP_orig partitioned from NEE_orig
compute_reco_gpp_from_nee_orig <- function(flux_data, ppfd_night_thr=10,
                                            gpp_init=c(0.08,15,0.5)) {
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
    ppfd_night_thr = ppfd_night_thr, gpp_init = gpp_init
  )
  flux_data$Reco_orig <- res$reco; flux_data$GPP_orig <- res$gpp
  log_msg("Reco_orig and GPP_orig computed from NEE_orig (ground truth, per regrowth period).")
  flux_data
}


# =============================================================================
# SECTION 7 — miniRECgap cross-validation — run_minirec_for_gap_size()
# =============================================================================
# The main cross-validation function iterates over every gap label within
# one gap-size category:
#
#   For each gap label:
#   1. Identify training rows: rows outside the gap with finite NEE.
#   2. Identify regrowth periods present in training rows.
#   3. For each regrowth period:
#      a. Fit R10 from nighttime training rows in that period.
#      b. Fit (α, GPPmax, θ) from daytime training rows in that period.
#      c. Apply the fitted parameters to ALL rows of that period to
#         compute Reco and GPP predictions.
#   4. Accumulate per-period Reco and GPP into full-length vectors.
#   5. Compute NEE_pred = Reco_pred - GPP_pred.
#   6. Write NEE_pred, Reco_pred, GPP_pred at the gap rows into the
#      respective prediction columns.
#
# Rows that belong to no named regrowth period (e.g. before the first
# grazing event of the year) are handled by a fallback global fit.
# Periods with fewer than 10 training rows receive NA predictions.

# PURPOSE
#   For each gap label, fit R10 and GPP parameters from outside-gap training
#   rows (per regrowth period), then for every gap row compute:
#     Reco_pred = R10 × A(T)
#     GPP_pred  = hyperbola(α, GPPmax, θ, PPFD)   [0 at night]
#     NEE_pred  = Reco_pred − GPP_pred
#
#   All three are stored at the gap rows:
#     NEE_{SIZE}_minirec_predicted   ← the primary CV target
#     Reco_{SIZE}_minirec_predicted  ← captured intermediate
#     GPP_{SIZE}_minirec_predicted   ← captured intermediate
#
#   This ensures NEE = Reco − GPP holds exactly and that Reco/GPP values
#   are the same mechanism that produced the NEE prediction — no separate
#   Reco/GPP model runs are needed.
#
# IN  : flux_data       — data frame with gap flags and masked NEE columns
#       gap_size_cat    — "S" | "M" | "L" | "VL"
#       ppfd_night_thr  — PPFD threshold: < = night, > = day  (default 10)

GPP_INIT_PARAMS <- c(alpha = 0.08, GPPmax = 15, theta = 0.5)

run_minirec_for_gap_size <- function(flux_data,
                                      gap_size_cat,
                                      ppfd_col       = "PPFD",
                                      temp_col       = "Temp",
                                      regrowth_col   = "regrowth_id",
                                      ppfd_night_thr = 10) {

  log_msg("miniRECgap [", gap_size_cat, "]: starting.")
  stopifnot(ppfd_col %in% names(flux_data),
            temp_col %in% names(flux_data),
            regrowth_col %in% names(flux_data))

  gap_labels <- names(flux_data)[grepl(paste0("^", gap_size_cat, "\\d+$"), names(flux_data))]
  gap_labels <- gap_labels[order(as.integer(sub(gap_size_cat, "", gap_labels)))]
  if (!length(gap_labels)) {
    log_msg("miniRECgap [", gap_size_cat, "]: no gap columns — skipping."); return(flux_data)
  }
  log_msg("miniRECgap [", gap_size_cat, "]: ", length(gap_labels), " gap labels.")

  nee_pred_col  <- paste0("NEE_",  gap_size_cat, "_minirec_predicted")
  reco_pred_col <- paste0("Reco_", gap_size_cat, "_minirec_predicted")
  gpp_pred_col  <- paste0("GPP_",  gap_size_cat, "_minirec_predicted")

  flux_data[[nee_pred_col]]  <- NA_real_
  flux_data[[reco_pred_col]] <- NA_real_
  flux_data[[gpp_pred_col]]  <- NA_real_

  for (gap_lbl in gap_labels) {
    masked_nee <- paste0("NEE_", gap_lbl)   # NA inside gap, observed outside
    if (!masked_nee %in% names(flux_data)) next

    ppfd     <- as.numeric(flux_data[[ppfd_col]])
    temp     <- as.numeric(flux_data[[temp_col]])
    nee_obs  <- as.numeric(flux_data[[masked_nee]])   # training signal
    regrowth <- flux_data[[regrowth_col]]
    gap_rows <- which(flux_data[[gap_lbl]] %in% c(TRUE, 1))
    is_obs   <- is.finite(nee_obs)
    is_night <- is.finite(ppfd) & ppfd < ppfd_night_thr
    is_day   <- is.finite(ppfd) & ppfd > ppfd_night_thr

    n <- nrow(flux_data)
    nee_pred_all  <- rep(NA_real_, n)
    reco_pred_all <- rep(NA_real_, n)
    gpp_pred_all  <- rep(NA_real_, n)

    regrowth_ids <- sort(unique(stats::na.omit(regrowth)))

    # ---- Fit per regrowth period ----------------------------------------
    for (rg_id in regrowth_ids) {
      rows_in_period  <- which(regrowth == rg_id)
      train_in_period <- intersect(rows_in_period, which(is_obs))
      if (length(train_in_period) < 10) next

      # Step 1 — R10 from nighttime training rows (NEE > 0)
      night_train <- train_in_period[is_night[train_in_period] &
                                       nee_obs[train_in_period] > 0]
      R10 <- if (length(night_train) > 0) {
        arr <- arrhenius_temperature_response(temp[night_train])
        denom <- sum(arr^2, na.rm = TRUE)
        if (is.finite(denom) && denom > 0)
          sum(nee_obs[night_train] * arr, na.rm = TRUE) / denom
        else 1.0
      } else 1.0

      # Reco for ALL rows in this period
      reco_period <- R10 * arrhenius_temperature_response(temp[rows_in_period])

      # Step 2 — GPP light-response from daytime training rows (NEE < 0)
      day_train <- train_in_period[is_day[train_in_period] &
                                     nee_obs[train_in_period] < 0]
      if (length(day_train) > 10) {
        reco_day_train <- R10 * arrhenius_temperature_response(temp[day_train])
        gpp_obs_train  <- reco_day_train - nee_obs[day_train]  # GPP = Reco − NEE

        opt <- tryCatch(
          optim(par    = GPP_INIT_PARAMS,
                fn     = gpp_residual_sum_squares,
                ppfd_vec = ppfd[day_train],
                gpp_obs  = gpp_obs_train,
                method   = "BFGS"),
          error = function(e) list(par = GPP_INIT_PARAMS)
        )
        gpp_params <- if (!is.null(opt$par)) opt$par else GPP_INIT_PARAMS
      } else {
        gpp_params <- GPP_INIT_PARAMS
      }

      # GPP for ALL rows in this period
      gpp_period <- light_response_gpp(gpp_params, ppfd[rows_in_period])
      gpp_period[!is_day[rows_in_period]]          <- 0        # GPP = 0 at night
      gpp_period[is.na(ppfd[rows_in_period])]      <- NA_real_

      # Accumulate predictions for all three variables
      reco_pred_all[rows_in_period] <- reco_period
      gpp_pred_all[rows_in_period]  <- gpp_period
      nee_pred_all[rows_in_period]  <- reco_period - gpp_period
    }

    # ---- Fallback: global fit for rows not covered by any period ----------
    need_fallback <- which(is.na(nee_pred_all) & is.finite(temp))
    if (length(need_fallback) > 0) {
      night_all <- which(is_night & is_obs & nee_obs > 0)
      if (length(night_all) >= 5) {
        arr_g <- arrhenius_temperature_response(temp[night_all])
        dg    <- sum(arr_g^2, na.rm = TRUE)
        R10g  <- if (is.finite(dg) && dg > 0)
          sum(nee_obs[night_all] * arr_g, na.rm = TRUE) / dg else 1.0
        reco_fb <- R10g * arrhenius_temperature_response(temp[need_fallback])

        day_all <- which(is_day & is_obs & nee_obs < 0)
        if (length(day_all) > 20) {
          reco_day_g <- R10g * arrhenius_temperature_response(temp[day_all])
          gpp_obs_g  <- reco_day_g - nee_obs[day_all]
          opt_g <- tryCatch(
            optim(par = GPP_INIT_PARAMS, fn = gpp_residual_sum_squares,
                  ppfd_vec = ppfd[day_all], gpp_obs = gpp_obs_g, method = "BFGS"),
            error = function(e) list(par = GPP_INIT_PARAMS)
          )
          gpp_p_g <- if (!is.null(opt_g$par)) opt_g$par else GPP_INIT_PARAMS
          gpp_fb  <- light_response_gpp(gpp_p_g, ppfd[need_fallback])
          gpp_fb[!is_day[need_fallback]]         <- 0
          gpp_fb[is.na(ppfd[need_fallback])]     <- NA_real_
          reco_pred_all[need_fallback] <- reco_fb
          gpp_pred_all[need_fallback]  <- gpp_fb
          nee_pred_all[need_fallback]  <- reco_fb - gpp_fb
        } else {
          # Day fit failed — only Reco available from fallback
          reco_pred_all[need_fallback] <- reco_fb
          nee_pred_all[need_fallback]  <- reco_fb   # GPP = 0 (no light data)
        }
      }
    }

    # ---- Write predictions to gap rows only --------------------------------
    flux_data[[nee_pred_col]][gap_rows]  <- nee_pred_all[gap_rows]
    flux_data[[reco_pred_col]][gap_rows] <- reco_pred_all[gap_rows]
    flux_data[[gpp_pred_col]][gap_rows]  <- gpp_pred_all[gap_rows]
  }

  log_msg("miniRECgap [", gap_size_cat, "] complete — ",
          "NEE MAE = ",  signif(calc_mae( flux_data$NEE_orig, flux_data[[nee_pred_col]]),3),
          "  RMSE = ",   signif(calc_rmse(flux_data$NEE_orig, flux_data[[nee_pred_col]]),3),
          "  R² = ",     signif(calc_r2(  flux_data$NEE_orig, flux_data[[nee_pred_col]]),3))
  flux_data
}


# =============================================================================
# SECTION 8 — Data preparation
# =============================================================================
# Filter to rows with finite NEE_orig.  Check availability of PPFD and Temp,
# which are required by both the Reco and GPP models.

req_cols <- c("timestamp", "NEE_orig", "PPFD", "Temp")
miss_cols <- setdiff(req_cols, names(df))
if (length(miss_cols)) stop("Missing required columns: ", paste(miss_cols, collapse = ", "))

if (!inherits(df$timestamp, "POSIXt"))
  df$timestamp <- suppressWarnings(
    lubridate::parse_date_time(df$timestamp,
      orders = c("Ymd HMS","Ymd HM","Ymd","Y/m/d HMS","Y/m/d"), tz = "UTC"))
stopifnot(!any(is.na(df$timestamp)))

flux_data <- df %>% dplyr::filter(is.finite(NEE_orig))
log_msg("Rows with finite NEE_orig retained: ", nrow(flux_data))


# =============================================================================
# SECTION 9 — Assign regrowth periods
# =============================================================================
# If Grazing_days_since is present in the data, regrowth periods are
# delineated by grazing events.  Otherwise a single global period is used.
# The regrowth_id column is passed to run_minirec_for_gap_size() to
# allow separate parameter estimation within each regrowth cycle.

if ("Grazing_days_since" %in% names(flux_data)) {
  flux_data <- flux_data %>%
    dplyr::mutate(regrowth_id = assign_regrowth_periods(Grazing_days_since, threshold = 1.0))
  log_msg("Regrowth periods: ", dplyr::n_distinct(flux_data$regrowth_id), " unique periods.")
} else {
  flux_data$regrowth_id <- 1L
  log_msg("Grazing_days_since absent — using one global regrowth period.")
}
# Compute ground-truth Reco_orig / GPP_orig from NEE_orig per regrowth period
if (all(c("PPFD","Temp") %in% names(flux_data)))
  flux_data <- compute_reco_gpp_from_nee_orig(flux_data)


# =============================================================================
# SECTION 10 — Construct artificial gaps
# =============================================================================
# Build the VL / L / M / S gap flag and masked-NEE matrices.
# See RF_CV.R Section 10 for gap-size specifications.

gap_cfg <- list(
  VL = list(base=30*48, tol=250, prefix="VL"),
  L  = list(base=14*48, tol=120, prefix="L"),
  M  = list(base= 7*48, tol= 60, prefix="M"),
  S  = list(base= 3*48, tol= 10, prefix="S")
)
log_msg("Constructing artificial gaps (VL → L → M → S) ...")
cv_VL <- create_cv_gaps(flux_data,   gap_cfg$VL$base, gap_cfg$VL$tol, "VL")
cv_L  <- create_cv_gaps(cv_VL$data,  gap_cfg$L$base,  gap_cfg$L$tol,  "L")
cv_M  <- create_cv_gaps(cv_L$data,   gap_cfg$M$base,  gap_cfg$M$tol,  "M")
cv_S  <- create_cv_gaps(cv_M$data,   gap_cfg$S$base,  gap_cfg$S$tol,  "S")
flux_data <- cv_S$data
log_msg("Gaps ready — VL:",length(cv_VL$sizes)," L:",length(cv_L$sizes),
        " M:",length(cv_M$sizes)," S:",length(cv_S$sizes))


# =============================================================================
# SECTION 11 — Run miniRECgap cross-validation  (all gap sizes)
# =============================================================================
# Apply run_minirec_for_gap_size() to all four categories.  Unlike RF/MLP/
# XGBoost, no separate post-hoc Reco/GPP derivation step is needed because
# miniRECgap produces Reco and GPP directly during gap prediction.

# Each call fills NEE, Reco, AND GPP prediction columns for that gap size.

log_msg("=== Running miniRECgap cross-validation ===")
for (cat in c("VL","L","M","S"))
  flux_data <- run_minirec_for_gap_size(flux_data, cat)
log_msg("All miniRECgap gap sizes complete.")


# =============================================================================
# SECTION 12 — Save outputs
# =============================================================================
# Save the augmented data frame and run metadata.

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
  paste("Targets         :", "NEE, Reco (intermediate), GPP (intermediate), Reco_orig, GPP_orig (ground truth)"),
  paste("Regrowth periods:", dplyr::n_distinct(flux_data$regrowth_id))
), file.path(RESULTS_DIR, "run_info.txt"))
log_msg("Saved: run_info.txt — miniRECgap script finished.")

# ========================== end miniRECgap_CV.R ================================
