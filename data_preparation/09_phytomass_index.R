# 09__phytomass_index.R — Phytomass Index (PI) per artificial gap
#
# For each site, loads the CV data frame from script 08 and appends one
# PI_{gap_label} column per gap label (across VL / L / M / S). PI is the
# normalised night-minus-day NEE contrast (Lohila 2004; Aurela 2001), computed
# only from observations OUTSIDE each gap to avoid leakage, then linearly
# interpolated across the gap rows. Thresholds: night PPFD < 1, day PPFD > 400;
# smoothing is a 21-day centred rolling pooled mean over a complete date spine.
#
# Input  : data/data_prepared/{SITE}_cv.rds   (from script 08)
# Output : data/data_prepared/{SITE}_cv.rds   (same file, PI columns appended)

library(here)
library(dplyr)
library(tibble)
library(tidyr)
library(zoo)

SITES   <- c("JC1", "JC2")
CV_DIR  <- here::here("data", "data_prepared")

NIGHT_PPFD  <- 1      # umol m-2 s-1 (Lohila 2004)
DAY_PPFD    <- 400    # umol m-2 s-1
WINDOW_DAYS <- 21L    # centred rolling window width (days)

# compute PI_{gap_label} for every gap label: build daily night/day sums from
# outside-gap rows, roll over an exact calendar spine, normalise to [0,1], then
# blank the gap rows and interpolate through them
compute_pi_per_gap <- function(df,
                               gap_prefixes = c("VL","L","M","S"),
                               night_ppfd   = NIGHT_PPFD,
                               day_ppfd     = DAY_PPFD,
                               window_days  = WINDOW_DAYS) {
  
  cd   <- as.Date(df$timestamp)
  nv   <- as.numeric(df$NEE_orig)
  pv   <- as.numeric(df$PPFD)
  ts_n <- as.numeric(df$timestamp)
  
  # complete calendar spine so the window always spans exact calendar days
  all_dates <- tibble::tibble(date = seq(min(cd), max(cd), by = "day"))
  
  rol <- function(x) zoo::rollapply(x, width = window_days, FUN = sum,
                                    align = "center", fill = NA_real_)
  
  for (pfx in gap_prefixes) {
    
    gap_cols <- names(df)[grepl(paste0("^", pfx, "\\d+$"), names(df))]
    gap_cols <- gap_cols[order(as.integer(sub(pfx, "", gap_cols)))]
    if (!length(gap_cols)) next
    
    for (gl in gap_cols) {
      
      # rows inside this gap vs outside
      gi      <- which(df[[gl]] %in% c(TRUE, 1))
      outside <- rep(TRUE, nrow(df)); outside[gi] <- FALSE
      
      ok_out  <- outside & is.finite(nv) & is.finite(pv)
      isn_out <- ok_out & pv < night_ppfd   # night half-hours outside gap
      isd_out <- ok_out & pv > day_ppfd     # day   half-hours outside gap
      
      # daily night/day sums and counts from outside-gap rows only
      dt <- tibble::tibble(
        date   = cd,
        nn     = dplyr::if_else(isn_out, nv, 0), cn = as.integer(isn_out),
        nd     = dplyr::if_else(isd_out, nv, 0), cd_col = as.integer(isd_out)
      ) |>
        dplyr::group_by(date) |>
        dplyr::summarise(nn = sum(nn), cn = sum(cn),
                         nd = sum(nd), cd = sum(cd_col), .groups = "drop") |>
        dplyr::arrange(date)
      
      # fill missing dates with zero so the window spans exact calendar days
      dt_full <- all_dates |>
        dplyr::left_join(dt, by = "date") |>
        dplyr::mutate(dplyr::across(c(nn, cn, nd, cd),
                                    ~ tidyr::replace_na(., 0)))
      
      # 21-day rolling pooled night mean (mn) and day mean (md)
      mn <- dplyr::if_else(rol(dt_full$cn) > 0,
                           rol(dt_full$nn) / rol(dt_full$cn), NA_real_)
      md <- dplyr::if_else(rol(dt_full$cd) > 0,
                           rol(dt_full$nd) / rol(dt_full$cd), NA_real_)
      PR <- dplyr::if_else(is.finite(mn) & is.finite(md), mn - md, NA_real_)
      mx <- suppressWarnings(max(PR, na.rm = TRUE))
      
      # normalise to [0, 1] and clamp
      PN <- if (!is.finite(mx) || mx <= 0) {
        rep(NA_real_, length(PR))
      } else {
        pmax(0, pmin(PR / mx, 1))
      }
      
      # map daily PI back to half-hourly rows
      ph <- PN[match(cd, dt_full$date)]
      
      # gap rows have no outside-gap observations — set to NA for interpolation
      ph[gi] <- NA_real_
      
      # linear interpolation through gap rows (rule = 2: extend flat at edges)
      ord <- order(df$timestamp)
      po  <- ph[ord]; xo <- ts_n[ord]
      rl  <- c(if (is.na(po[1])) 2 else 1,
               if (is.na(po[length(po)])) 2 else 1)
      if (sum(is.finite(po)) >= 2)
        po <- zoo::na.approx(po, x = xo, na.rm = FALSE, rule = rl)
      ph[ord] <- po
      
      df[[paste0("PI_", gl)]] <- ph
    }
  }
  
  df
}

for (site in SITES) {
  
  # load the CV data frame from script 08 (skip the site if missing)
  cv_path <- file.path(CV_DIR, paste0(site, "_cv.rds"))
  if (!file.exists(cv_path)) {
    warning("CV file not found, skipping: ", cv_path); next
  }
  df <- readRDS(cv_path)
  
  # required columns for PI
  req <- c("timestamp", "NEE_orig", "PPFD")
  missing_req <- setdiff(req, names(df))
  if (length(missing_req))
    stop("Required columns missing from ", cv_path, ": ",
         paste(missing_req, collapse = ", "))
  
  # compute and append the PI columns
  df <- compute_pi_per_gap(df)
  
  # overwrite the CV file with PI columns added
  saveRDS(df, cv_path)
}
