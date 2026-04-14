# =============================================================================
# SCRIPT 09 — Phytomass Index (PI) per Artificial Gap
# =============================================================================
# PURPOSE
#   For each site (JC1, JC2), load the CV data frame produced by Script 08
#   (JCi_cv.rds) and compute one PI_{gap_label} column per gap label across
#   all four gap-size categories (VL / L / M / S).
#
#   PI is an empirically determined proxy for above-ground phytomass, defined
#   as the normalised contrast between nighttime and daytime mean NEE
#   (Lohila et al. 2004; Aurela et al. 2001):
#
#     PR(d) = mean_night_NEE(d) − mean_day_NEE(d)
#     PI(d) = PR(d) / max_d( PR(d) )    ∈ [0, 1]
#
#   To prevent information leakage, PI for gap label Si is computed using
#   only observations OUTSIDE Si.  Gap rows (Si = TRUE) are set to NA and
#   then recovered by linear interpolation between the nearest valid PI
#   values on either side of the gap.
#
#   PPFD thresholds (Lohila 2004):
#     Night : PPFD < 1   µmol m⁻² s⁻¹
#     Day   : PPFD > 400 µmol m⁻² s⁻¹
#
#   Smoothing: 21-day centred rolling pooled mean.  A complete calendar date
#   spine is used so the window always spans exactly 21 calendar days,
#   regardless of data gaps.
#
# INPUTS
#   data/data_prepared/JC1_cv.rds   — output of Script 08
#   data/data_prepared/JC2_cv.rds
#
# OUTPUTS
#   data/data_prepared/JC1_cv.rds   — same file, overwritten with PI columns
#   data/data_prepared/JC2_cv.rds
#
#   New columns added (one per gap label):
#     PI_VL1, PI_VL2, …   — PI computed outside each VL gap
#     PI_L1,  PI_L2,  …
#     PI_M1,  PI_M2,  …
#     PI_S1,  PI_S2,  …
# =============================================================================

library(here)
library(dplyr)
library(tibble)
library(tidyr)
library(zoo)

# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------
SITES   <- c("JC1", "JC2")
CV_DIR  <- here::here("data", "data_prepared")

NIGHT_PPFD  <- 1      # µmol m⁻² s⁻¹  (Lohila 2004)
DAY_PPFD    <- 400    # µmol m⁻² s⁻¹
WINDOW_DAYS <- 21L    # centred rolling window width (days)

# -----------------------------------------------------------------------------
# HELPER — compute PI_{gap_label} for every gap label in df
# -----------------------------------------------------------------------------
# For each gap label gl (e.g. "M3"):
#   1. Identify gap rows (gl == TRUE).
#   2. Build daily sums/counts from OUTSIDE-gap rows only.
#   3. Join to a complete calendar spine so the rolling window always spans
#      exactly WINDOW_DAYS calendar days (missing dates → sum/count = 0).
#   4. Compute 21-day rolling pooled night mean (mn) and day mean (md).
#   5. PR = mn − md;  normalise to [0,1] → PI.
#   6. Set gap rows to NA, then interpolate linearly using zoo::na.approx.

compute_pi_per_gap <- function(df,
                               gap_prefixes = c("VL","L","M","S"),
                               night_ppfd   = NIGHT_PPFD,
                               day_ppfd     = DAY_PPFD,
                               window_days  = WINDOW_DAYS) {
  
  cd   <- as.Date(df$timestamp)
  nv   <- as.numeric(df$NEE_orig)
  pv   <- as.numeric(df$PPFD)
  ts_n <- as.numeric(df$timestamp)
  
  # Complete calendar spine — rolling window always spans exact calendar days
  all_dates <- tibble::tibble(date = seq(min(cd), max(cd), by = "day"))
  
  rol <- function(x) zoo::rollapply(x, width = window_days, FUN = sum,
                                    align = "center", fill = NA_real_)
  
  for (pfx in gap_prefixes) {
    
    gap_cols <- names(df)[grepl(paste0("^", pfx, "\\d+$"), names(df))]
    gap_cols <- gap_cols[order(as.integer(sub(pfx, "", gap_cols)))]
    if (!length(gap_cols)) next
    
    message("  Computing PI for ", length(gap_cols),
            " '", pfx, "' gap labels ...")
    
    for (gl in gap_cols) {
      
      # Row indices inside this gap
      gi      <- which(df[[gl]] %in% c(TRUE, 1))
      outside <- rep(TRUE, nrow(df)); outside[gi] <- FALSE
      
      ok_out  <- outside & is.finite(nv) & is.finite(pv)
      isn_out <- ok_out & pv < night_ppfd   # night half-hours outside gap
      isd_out <- ok_out & pv > day_ppfd     # day   half-hours outside gap
      
      # Daily sums and counts from outside-gap rows only
      dt <- tibble::tibble(
        date   = cd,
        nn     = dplyr::if_else(isn_out, nv, 0), cn = as.integer(isn_out),
        nd     = dplyr::if_else(isd_out, nv, 0), cd_col = as.integer(isd_out)
      ) |>
        dplyr::group_by(date) |>
        dplyr::summarise(nn = sum(nn), cn = sum(cn),
                         nd = sum(nd), cd = sum(cd_col), .groups = "drop") |>
        dplyr::arrange(date)
      
      # Fill missing dates with zero so the window spans exact calendar days
      dt_full <- all_dates |>
        dplyr::left_join(dt, by = "date") |>
        dplyr::mutate(dplyr::across(c(nn, cn, nd, cd),
                                    ~ tidyr::replace_na(., 0)))
      
      # 21-day rolling pooled means
      mn <- dplyr::if_else(rol(dt_full$cn) > 0,
                           rol(dt_full$nn) / rol(dt_full$cn), NA_real_)
      md <- dplyr::if_else(rol(dt_full$cd) > 0,
                           rol(dt_full$nd) / rol(dt_full$cd), NA_real_)
      PR <- dplyr::if_else(is.finite(mn) & is.finite(md), mn - md, NA_real_)
      mx <- suppressWarnings(max(PR, na.rm = TRUE))
      
      # Normalise to [0, 1] and clamp
      PN <- if (!is.finite(mx) || mx <= 0) {
        rep(NA_real_, length(PR))
      } else {
        pmax(0, pmin(PR / mx, 1))
      }
      
      # Map daily PI back to half-hourly rows
      ph <- PN[match(cd, dt_full$date)]
      
      # Gap rows have no outside-gap observations — set to NA for interpolation
      ph[gi] <- NA_real_
      
      # Linear interpolation through gap rows (rule = 2: extend flat at edges)
      ord <- order(df$timestamp)
      po  <- ph[ord]; xo <- ts_n[ord]
      rl  <- c(if (is.na(po[1])) 2 else 1,
               if (is.na(po[length(po)])) 2 else 1)
      if (sum(is.finite(po)) >= 2)
        po <- zoo::na.approx(po, x = xo, na.rm = FALSE, rule = rl)
      ph[ord] <- po
      
      df[[paste0("PI_", gl)]] <- ph
      
      message("    PI_", gl,
              " | gap rows: ", length(gi),
              " | PI range: [", round(min(ph, na.rm = TRUE), 3),
              ", ",             round(max(ph, na.rm = TRUE), 3), "]",
              " | NA remaining: ", sum(is.na(ph)))
    }
  }
  
  df
}

# -----------------------------------------------------------------------------
# MAIN LOOP — process each site
# -----------------------------------------------------------------------------
for (site in SITES) {
  
  message("\n", strrep("=", 70))
  message("Processing site: ", site)
  message(strrep("=", 70))
  
  # --- 1. Load CV data frame from Script 08 ----------------------------------
  cv_path <- file.path(CV_DIR, paste0(site, "_cv.rds"))
  if (!file.exists(cv_path)) {
    warning("CV file not found, skipping: ", cv_path); next
  }
  df <- readRDS(cv_path)
  message("Loaded: ", cv_path, "  (", nrow(df), " rows × ", ncol(df), " cols)")
  
  # Verify required columns
  req <- c("timestamp", "NEE_orig", "PPFD")
  missing_req <- setdiff(req, names(df))
  if (length(missing_req))
    stop("Required columns missing from ", cv_path, ": ",
         paste(missing_req, collapse = ", "))
  
  # --- 2. Compute PI per gap label -------------------------------------------
  message("Computing PI columns ...")
  df <- compute_pi_per_gap(df)
  
  # --- 3. Report ---------------------------------------------------------------
  pi_cols <- names(df)[grepl("^PI_", names(df))]
  message("PI columns added: ", length(pi_cols))
  message("Output dimensions: ", nrow(df), " rows × ", ncol(df), " cols")
  
  # --- 4. Overwrite CV file with PI columns added ----------------------------
  saveRDS(df, cv_path)
  message("Saved (updated): ", cv_path)
}

message("\nScript 09 complete.")