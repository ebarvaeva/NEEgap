# =============================================================================
# miniRECgap — Cross-Validation Gap-Filling for NEE, Reco, and GPP
# =============================================================================
#
# WHAT THIS SCRIPT DOES
#   Estimates NEE inside each artificial gap using the miniRECgap process model:
#     NEE_pred = Reco_model − GPP_model
#   where Reco and GPP are parametric functions fitted from outside-gap
#   training observations per vegetation regrowth period.
#
# HOW RECO AND GPP ARE STORED
#   During run_minirec_for_gap_size(), Reco and GPP are computed as
#   intermediates alongside NEE_pred and stored directly:
#     Reco_{SIZE}_minirec_predicted  — captured intermediate
#     GPP_{SIZE}_minirec_predicted   — captured intermediate
#   This ensures NEE = Reco − GPP holds exactly at every gap row.
#
# INPUTS
#   df — a pre-built JC1_cv.rds / JC2_cv.rds data frame, sourced externally.
#        Already filtered to rows with finite NEE_orig and augmented with
#        gap flag columns (VL1…, L1…, M1…, S1…) and masked NEE columns.
#
# OUTPUTS  (written to RESULTS_DIR)
#   df_cv_all_predictions.rds
#       NEE_{S|M|L|VL}_minirec_predicted
#       Reco_{S|M|L|VL}_minirec_predicted
#       GPP_{S|M|L|VL}_minirec_predicted
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


# =============================================================================
# SECTION 3 — Progress logger
# =============================================================================
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
log_msg("miniRECgap started.  RESULTS_DIR = ", RESULTS_DIR)


# =============================================================================
# SECTION 4 — Data validation
# =============================================================================
req_cols <- c("timestamp", "NEE_orig", "PPFD", "Temp")
miss_cols <- setdiff(req_cols, names(df))
if (length(miss_cols)) stop("Missing required columns: ", paste(miss_cols, collapse = ", "))

if (!inherits(df$timestamp, "POSIXt"))
  df$timestamp <- suppressWarnings(
    lubridate::parse_date_time(df$timestamp,
                               orders = c("Ymd HMS", "Ymd HM", "Ymd", "Y/m/d HMS", "Y/m/d"), tz = "UTC"))
stopifnot(!any(is.na(df$timestamp)))

df <- dplyr::arrange(df, timestamp)
log_msg("Input rows: ", nrow(df), "  (pre-filtered to finite NEE_orig by Script 08)")


# =============================================================================
# SECTION 5 — Regrowth period detection
# =============================================================================
assign_regrowth_periods <- function(days_since, threshold = 1.0) {
  x       <- suppressWarnings(as.numeric(days_since))
  is_zero <- is.finite(x) & (x == 0); x_prev <- dplyr::lag(x)
  cumsum(ifelse(is_zero & (is.na(x_prev) | !is.finite(x_prev) | (x_prev > threshold)), 1L, 0L))
}

if ("Grazing_days_since" %in% names(df)) {
  df <- df %>%
    dplyr::mutate(regrowth_id = assign_regrowth_periods(Grazing_days_since, threshold = 1.0))
  log_msg("Regrowth periods: ", dplyr::n_distinct(df$regrowth_id), " unique periods.")
} else {
  df$regrowth_id <- 1L
  log_msg("Grazing_days_since absent — using one global regrowth period.")
}


# =============================================================================
# SECTION 6 — Process-based flux models
# =============================================================================
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

gpp_residual_sum_squares <- function(params, ppfd_vec, gpp_obs)
  sum((light_response_gpp(params, ppfd_vec) - gpp_obs)^2)

# Ground truth: Reco_orig and GPP_orig partitioned from NEE_orig
compute_reco_gpp_from_nee_orig <- function(flux_data, ppfd_night_thr = 10,
                                           gpp_init = c(0.08, 15, 0.5)) {
  if (!all(c("PPFD", "Temp", "regrowth_id", "NEE_orig") %in% names(flux_data))) {
    log_msg("compute_reco_gpp_from_nee_orig: required columns missing — skipping.")
    return(flux_data)
  }
  n        <- nrow(flux_data)
  ppfd_vec <- as.numeric(flux_data$PPFD)
  temp_vec <- as.numeric(flux_data$Temp)
  nee_vec  <- as.numeric(flux_data$NEE_orig)
  regrowth <- flux_data$regrowth_id
  is_night <- is.finite(ppfd_vec) & ppfd_vec < ppfd_night_thr
  is_day   <- is.finite(ppfd_vec) & ppfd_vec > ppfd_night_thr
  reco_out <- rep(NA_real_, n)
  gpp_out  <- rep(NA_real_, n)
  
  for (rg_id in sort(unique(stats::na.omit(regrowth)))) {
    rip        <- which(regrowth == rg_id)
    train_rows <- rip[is.finite(nee_vec[rip])]
    if (length(train_rows) < 10) next
    
    night_rows <- train_rows[is_night[train_rows] & nee_vec[train_rows] > 0]
    R10 <- if (length(night_rows) > 0) {
      arr <- arrhenius_temperature_response(temp_vec[night_rows])
      d   <- sum(arr^2, na.rm = TRUE)
      if (is.finite(d) && d > 0) sum(nee_vec[night_rows] * arr, na.rm = TRUE) / d else 1.0
    } else 1.0
    
    reco_out[rip] <- R10 * arrhenius_temperature_response(temp_vec[rip])
    
    day_rows <- train_rows[is_day[train_rows] & nee_vec[train_rows] < 0]
    if (length(day_rows) > 10) {
      gpp_obs <- R10 * arrhenius_temperature_response(temp_vec[day_rows]) - nee_vec[day_rows]
      opt <- tryCatch(
        optim(par = gpp_init, fn = gpp_residual_sum_squares,
              ppfd_vec = ppfd_vec[day_rows], gpp_obs = gpp_obs, method = "BFGS"),
        error = function(e) list(par = gpp_init))
      gpp_params <- if (!is.null(opt$par)) opt$par else gpp_init
    } else gpp_params <- gpp_init
    
    gpp_period        <- light_response_gpp(gpp_params, ppfd_vec[rip])
    gpp_period[!is_day[rip]]         <- 0
    gpp_period[is.na(ppfd_vec[rip])] <- NA_real_
    gpp_out[rip] <- gpp_period
  }
  
  flux_data$Reco_orig <- reco_out
  flux_data$GPP_orig  <- gpp_out
  log_msg("Reco_orig and GPP_orig computed from NEE_orig (ground truth, per regrowth period).")
  flux_data
}

df <- compute_reco_gpp_from_nee_orig(df)


# =============================================================================
# SECTION 7 — miniRECgap cross-validation
# =============================================================================
# For each gap label, fit R10 and GPP parameters from outside-gap training
# rows (per regrowth period), then for every gap row compute:
#   Reco_pred = R10 × A(T)
#   GPP_pred  = hyperbola(α, GPPmax, θ, PPFD)   [0 at night]
#   NEE_pred  = Reco_pred − GPP_pred
#
# All three are written to gap rows only:
#   NEE_{SIZE}_minirec_predicted
#   Reco_{SIZE}_minirec_predicted
#   GPP_{SIZE}_minirec_predicted

GPP_INIT_PARAMS <- c(alpha = 0.08, GPPmax = 15, theta = 0.5)

run_minirec_for_gap_size <- function(flux_data, gap_size_cat,
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
    masked_nee <- paste0("NEE_", gap_lbl)
    if (!masked_nee %in% names(flux_data)) next
    
    ppfd     <- as.numeric(flux_data[[ppfd_col]])
    temp     <- as.numeric(flux_data[[temp_col]])
    nee_obs  <- as.numeric(flux_data[[masked_nee]])
    regrowth <- flux_data[[regrowth_col]]
    gap_rows <- which(flux_data[[gap_lbl]] %in% c(TRUE, 1))
    is_obs   <- is.finite(nee_obs)
    is_night <- is.finite(ppfd) & ppfd < ppfd_night_thr
    is_day   <- is.finite(ppfd) & ppfd > ppfd_night_thr
    
    n             <- nrow(flux_data)
    nee_pred_all  <- rep(NA_real_, n)
    reco_pred_all <- rep(NA_real_, n)
    gpp_pred_all  <- rep(NA_real_, n)
    
    for (rg_id in sort(unique(stats::na.omit(regrowth)))) {
      rows_in_period  <- which(regrowth == rg_id)
      train_in_period <- intersect(rows_in_period, which(is_obs))
      if (length(train_in_period) < 10) next
      
      # Step 1 — R10 from nighttime training rows (NEE > 0)
      night_train <- train_in_period[is_night[train_in_period] & nee_obs[train_in_period] > 0]
      R10 <- if (length(night_train) > 0) {
        arr   <- arrhenius_temperature_response(temp[night_train])
        denom <- sum(arr^2, na.rm = TRUE)
        if (is.finite(denom) && denom > 0)
          sum(nee_obs[night_train] * arr, na.rm = TRUE) / denom
        else 1.0
      } else 1.0
      
      reco_period <- R10 * arrhenius_temperature_response(temp[rows_in_period])
      
      # Step 2 — GPP light-response from daytime training rows (NEE < 0)
      day_train <- train_in_period[is_day[train_in_period] & nee_obs[train_in_period] < 0]
      if (length(day_train) > 10) {
        gpp_obs_train <- R10 * arrhenius_temperature_response(temp[day_train]) - nee_obs[day_train]
        opt <- tryCatch(
          optim(par = GPP_INIT_PARAMS, fn = gpp_residual_sum_squares,
                ppfd_vec = ppfd[day_train], gpp_obs = gpp_obs_train, method = "BFGS"),
          error = function(e) list(par = GPP_INIT_PARAMS))
        gpp_params <- if (!is.null(opt$par)) opt$par else GPP_INIT_PARAMS
      } else gpp_params <- GPP_INIT_PARAMS
      
      gpp_period <- light_response_gpp(gpp_params, ppfd[rows_in_period])
      gpp_period[!is_day[rows_in_period]]         <- 0
      gpp_period[is.na(ppfd[rows_in_period])]     <- NA_real_
      
      reco_pred_all[rows_in_period] <- reco_period
      gpp_pred_all[rows_in_period]  <- gpp_period
      nee_pred_all[rows_in_period]  <- reco_period - gpp_period
    }
    
    # Fallback: global fit for rows not covered by any regrowth period
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
          gpp_obs_g <- R10g * arrhenius_temperature_response(temp[day_all]) - nee_obs[day_all]
          opt_g <- tryCatch(
            optim(par = GPP_INIT_PARAMS, fn = gpp_residual_sum_squares,
                  ppfd_vec = ppfd[day_all], gpp_obs = gpp_obs_g, method = "BFGS"),
            error = function(e) list(par = GPP_INIT_PARAMS))
          gpp_p_g <- if (!is.null(opt_g$par)) opt_g$par else GPP_INIT_PARAMS
          gpp_fb  <- light_response_gpp(gpp_p_g, ppfd[need_fallback])
          gpp_fb[!is_day[need_fallback]]     <- 0
          gpp_fb[is.na(ppfd[need_fallback])] <- NA_real_
          reco_pred_all[need_fallback] <- reco_fb
          gpp_pred_all[need_fallback]  <- gpp_fb
          nee_pred_all[need_fallback]  <- reco_fb - gpp_fb
        } else {
          reco_pred_all[need_fallback] <- reco_fb
          nee_pred_all[need_fallback]  <- reco_fb
        }
      }
    }
    
    flux_data[[nee_pred_col]][gap_rows]  <- nee_pred_all[gap_rows]
    flux_data[[reco_pred_col]][gap_rows] <- reco_pred_all[gap_rows]
    flux_data[[gpp_pred_col]][gap_rows]  <- gpp_pred_all[gap_rows]
    log_msg("miniRECgap [", gap_size_cat, "] ", gap_lbl, ": filled ", length(gap_rows), " gap rows.")
  }
  
  log_msg("miniRECgap [", gap_size_cat, "] complete.")
  flux_data
}


# =============================================================================
# SECTION 8 — Run miniRECgap for all gap sizes
# =============================================================================
log_msg("=== Running miniRECgap cross-validation ===")
for (cat in c("VL", "L", "M", "S"))
  df <- run_minirec_for_gap_size(df, cat)
log_msg("All miniRECgap gap sizes complete.")


# =============================================================================
# SECTION 9 — Save output
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
  paste("Targets         :", "NEE, Reco (intermediate), GPP (intermediate), Reco_orig, GPP_orig (ground truth)"),
  paste("Regrowth periods:", dplyr::n_distinct(df$regrowth_id))
), file.path(RESULTS_DIR, "run_info.txt"))
log_msg("Saved: run_info.txt — miniRECgap script finished.")
# ========================== end miniRECgap_CV.R ================================