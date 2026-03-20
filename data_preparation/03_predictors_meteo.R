# =============================================================================
# 03_predictors_meteo.R — Meteorological Predictors and Temporal Encodings
# =============================================================================
#
# PURPOSE
#   Builds the full BASE predictor set from the full-year QC grid (Script 02)
#   and external meteorological data sources.  The following columns are
#   created or gap-filled:
#
#   PPFD  -- photosynthetic photon flux density (umol m-2 s-1)
#     Priority 1: on-site biomet (PPFD_1_1_1_biomet or PPFD_1_1_1).
#     Priority 2: twin-site biomet (PPFD from the other grassland site).
#     Priority 3: Rg * 0.45 * 4.57 (PAR fraction x photon energy factor;
#                 McCree 1972).  Applied only where PPFD is still NA after
#                 Priorities 1 and 2.
#     Abnormal PPFD periods (< 0 or >= 250 umol m-2 s-1 for >= 3 consecutive
#     days, or all < 375 umol m-2 s-1 for >= 7 days) are set to NA before
#     fallback imputation to avoid propagating sensor malfunctions.
#
#   Rg    -- global shortwave radiation (W m-2), Met Eirean solar hourly,
#            linearly interpolated from hourly to 30-min.
#
#   Temp  -- air temperature (degC), Met Eirean hourly, linearly interpolated.
#   RH    -- relative humidity (%), Met Eirean hourly, linearly interpolated.
#   rain  -- precipitation (mm per half-hour), Met Eirean hourly, interpolated.
#   VPD   -- vapour pressure deficit (kPa), derived from Temp and RH via
#            the Tetens equation (Magnus approximation).
#
#   rain_rolling_24  -- 24-h trailing mean of rain (right-aligned, 48 rows).
#                       Captures antecedent soil moisture conditions.
#   night            -- binary flag: 1 = PPFD < 10 umol m-2 s-1, else 0.
#
#   Temporal encodings (remove circular discontinuity at period boundaries):
#     hour_sin, hour_cos  -- hour of day encoded as sine/cosine pair
#     doy_sin,  doy_cos   -- day of year encoded as sine/cosine pair
#     month_sin,month_cos -- month encoded as sine/cosine pair
#   Season dummies: Winter (Nov-Jan), Spring (Feb-Apr), Summer (May-Jul),
#                   Autumn (Aug-Oct) -- Irish meteorological calendar.
#
#   NEE quality filter (applied last):
#     Night-time NEE < 0 is physiologically implausible (plants cannot
#     photosynthesize without light) and is set to NA.  The u* filter
#     was applied in Script 01.
#
# INPUTS
#   data/data_qc/{SITE}_{YEAR}/merged_qc_fullgrid.csv  (Script 02 output)
#   data/meteireann_data/met_hourly.csv                 hourly Temp, RH, rain
#   data/meteireann_data/solar_hourly.xlsx              hourly global radiation
#   data/raw_data/{TWIN_SITE}_{YEAR}_biomet.xlsx        twin-site PPFD fallback
#
# OUTPUTS
#   data/data_prepared/{SITE}_{YEAR}_predictors.rds
#
# SECTIONS
#   User settings
#   Helper  detect_bad_periods() / apply_bad_periods_na() -- PPFD QC
#   1  Load full-grid QC data     -- standardise NEE and PPFD column names
#   2  Twin-site PPFD fallback    -- fill on-site PPFD gaps from other site
#   3  Rg-derived PPFD fallback   -- described in Purpose above
#   4  Coalesce PPFD              -- merge Priority 1 and 2; Priority 3 in S6
#   5  Met Eirean: Temp, RH, rain -- load and interpolate to 30 min; derive VPD
#   6  Met Eirean: Rg             -- load solar, interpolate, apply Rg->PPFD fill
#   7  Temporal encodings         -- sine/cosine pairs and season dummies
#   8  Rolling rain and night     -- rain_rolling_24, night flag, NEE filter
#   9  Select and save
#
# =============================================================================

library(here)
library(dplyr)
library(readr)
library(readxl)
library(lubridate)
library(tidyr)
library(zoo)
library(ranger)

# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------
# TWIN_SITE is set automatically to the other site; change if needed.
if (!exists("TWIN_SITE")) TWIN_SITE <- setdiff(c("JC1", "JC2"), SITE)
TZ        <- "UTC"

# PPFD abnormal period detection thresholds.
# A 'bad period' is flagged when PPFD is outside [0, 250) for >= PPFD_BADPERIOD_DAYS
# consecutive days, or when all values are below PPFD_ALLLT_THRESHOLD for
# >= PPFD_ALLLT_DAYS days (indicating a covered sensor or calibration failure).
# These values can be adjusted if sensor behaviour at a new site differs.
PPFD_BADPERIOD_DAYS    <- 3     # min consecutive days outside [0, 250) to flag
PPFD_ALLLT_THRESHOLD   <- 375   # umol m-2 s-1 -- 'all below' threshold
PPFD_ALLLT_DAYS        <- 7     # min consecutive days all-below threshold

OUT_DIR <- here::here("data", "data_prepared")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# HELPER: detect_bad_periods() and apply_bad_periods_na()
# -----------------------------------------------------------------------------
# detect_bad_periods() identifies contiguous runs of suspicious PPFD values
# using run-length encoding (rle()).  Two rules are applied independently:
#   Rule A: PPFD outside [0, 250) for >= days * interval_per_day half-hours.
#   Rule B: All PPFD < all_lt_threshold for >= all_lt_days * interval_per_day.
# Detected segments are merged if they overlap or are adjacent.
# apply_bad_periods_na() sets the flagged half-hours to NA so that the
# fallback chain (twin site, Rg formula) fills them rather than propagating
# the bad sensor readings.
detect_bad_periods <- function(df, col,
                                days           = 3L,
                                interval_per_day = 48L,
                                all_lt_threshold = 375,
                                all_lt_days      = 7L) {
  stopifnot(col %in% names(df))
  df <- dplyr::arrange(df, timestamp)
  x  <- df[[col]]

  get_segments <- function(flag, min_len) {
    if (!length(flag))
      return(tibble::tibble(start_idx = integer(), end_idx = integer()))
    r      <- rle(flag)
    ends   <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1L
    keep   <- which(r$values & r$lengths >= min_len)
    if (!length(keep))
      return(tibble::tibble(start_idx = integer(), end_idx = integer()))
    tibble::tibble(start_idx = starts[keep], end_idx = ends[keep])
  }

  segA <- get_segments(!is.na(x) & (x < 0 | x >= 250 | is.na(x)),
                        as.integer(days * interval_per_day))
  segB <- get_segments(!is.na(x) & x < all_lt_threshold,
                        as.integer(all_lt_days * interval_per_day))

  segs <- dplyr::bind_rows(segA, segB) %>% dplyr::arrange(start_idx, end_idx)
  if (!nrow(segs)) return(tibble::tibble(start_idx = integer(), end_idx = integer(),
                                          start_time = as.POSIXct(character()),
                                          end_time   = as.POSIXct(character())))

  merged <- list(); cur <- segs[1L, ]
  for (i in seq_len(nrow(segs))[-1]) {
    nxt <- segs[i, ]
    if (nxt$start_idx <= cur$end_idx + 1L)
      cur$end_idx <- max(cur$end_idx, nxt$end_idx)
    else { merged[[length(merged)+1L]] <- cur; cur <- nxt }
  }
  merged[[length(merged)+1L]] <- cur
  dplyr::bind_rows(merged) %>%
    dplyr::mutate(start_time = df$timestamp[start_idx],
                  end_time   = df$timestamp[end_idx])
}

apply_bad_periods_na <- function(df, col, segs) {
  if (!nrow(segs)) return(df)
  bad <- rep(FALSE, nrow(df))
  for (i in seq_len(nrow(segs)))
    bad[segs$start_idx[i]:segs$end_idx[i]] <- TRUE
  df[[col]][bad] <- NA_real_
  df
}

# -----------------------------------------------------------------------------
# SECTION 1 — Load full-grid QC data
# -----------------------------------------------------------------------------
# Reads the output of Script 02 and standardises two column names that vary
# across EddyPro versions and QC workflows:
#   co2_flux_base_filters_6_70_grid -> renamed to NEE_orig
#   PPFD_1_1_1_biomet or PPFD_1_1_1 -> renamed to PPFD_biomet
# Bad PPFD periods are set to NA using the helper defined above so that
# the fallback imputation in Sections 2-4 does not perpetuate sensor errors.
message("Loading QC full-grid data ...")
in_path <- here::here("data", "data_qc",
                       paste0("results_qc_", SITE, "_", YEAR), "merged_qc_fullgrid.csv")
if (!file.exists(in_path)) stop("Input not found: ", in_path)

df <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(timestamp = as.POSIXct(timestamp, tz = TZ)) %>%
  arrange(timestamp) %>%
  distinct(timestamp, .keep_all = TRUE)

# Standardise NEE column name
if (!"NEE_orig" %in% names(df) && "co2_flux_base_filters_6_70_grid" %in% names(df))
  df <- df %>% rename(NEE_orig = co2_flux_base_filters_6_70_grid)

# Standardise on-site PPFD column
ppfd_biomet_col <- intersect(c("PPFD_1_1_1_biomet", "PPFD_1_1_1"), names(df))[1]
if (is.na(ppfd_biomet_col)) {
  df$PPFD_biomet <- NA_real_
} else {
  df <- df %>%
    rename(PPFD_biomet = all_of(ppfd_biomet_col)) %>%
    mutate(PPFD_biomet = na_if(PPFD_biomet, -9999),
           PPFD_biomet = if_else(!is.na(PPFD_biomet) & PPFD_biomet < 0,
                                  -PPFD_biomet, PPFD_biomet))
}

# Apply bad-period filter to on-site PPFD
bp <- detect_bad_periods(df, "PPFD_biomet",
                          days           = PPFD_BADPERIOD_DAYS,
                          all_lt_threshold = PPFD_ALLLT_THRESHOLD,
                          all_lt_days    = PPFD_ALLLT_DAYS)
df <- apply_bad_periods_na(df, "PPFD_biomet", bp)
message("  On-site PPFD NA after bad-period filter: ", sum(is.na(df$PPFD_biomet)))

# -----------------------------------------------------------------------------
# SECTION 2 — Twin-site PPFD fallback
# -----------------------------------------------------------------------------
# When on-site PPFD is unavailable (NA after bad-period filtering), the biomet
# file from the other grassland site (TWIN_SITE) is used as a surrogate.
# The two sites are approximately 2 km apart, so their PPFD is highly correlated
# on the half-hourly timescale.  The twin-site PPFD is also subjected to the
# same bad-period filter before use to prevent propagating its own sensor errors.
# The left_join on timestamp aligns the two timeseries exactly.
message("Loading twin-site PPFD fallback from biomet (", TWIN_SITE, ") ...")

twin_ppfd_col <- paste0("PPFD_alt_", TWIN_SITE)
twin_biomet_path <- here::here("data", "raw_data",
                                paste0(TWIN_SITE, "_", YEAR, "_biomet.xlsx"))

if (file.exists(twin_biomet_path)) {
  twin_df <- read_excel(twin_biomet_path)%>%
    distinct(date, time, .keep_all = TRUE) %>%
    mutate(
      timestamp = update(ymd_hms(time, tz = "Europe/Dublin"),
                         year  = year(as_date(date)),
                         month = month(as_date(date)),
                         mday  = mday(as_date(date)))
    ) %>%
    arrange(timestamp) %>%
    distinct(timestamp, .keep_all = TRUE) %>%
    select(timestamp, any_of(c("PPFD_1_1_1_biomet", "PPFD_1_1_1"))) %>%
    rename_with(~ twin_ppfd_col, -timestamp) %>%
    mutate(
      !!twin_ppfd_col := na_if(.data[[twin_ppfd_col]], -9999),
      !!twin_ppfd_col := if_else(
        !is.na(.data[[twin_ppfd_col]]) & .data[[twin_ppfd_col]] < 0,
        -.data[[twin_ppfd_col]], .data[[twin_ppfd_col]])
    )

  df <- df %>% left_join(twin_df, by = "timestamp")
  bp2 <- detect_bad_periods(df, twin_ppfd_col,
                              days           = PPFD_BADPERIOD_DAYS,
                              all_lt_threshold = PPFD_ALLLT_THRESHOLD,
                              all_lt_days    = PPFD_ALLLT_DAYS)
  df <- apply_bad_periods_na(df, twin_ppfd_col, bp2)
  message("  Twin-site PPFD loaded: ", sum(is.finite(df[[twin_ppfd_col]])),
          " finite values")
} else {
  message("  Twin-site biomet not found: ", twin_biomet_path, " — skipping.")
  df[[twin_ppfd_col]] <- NA_real_
}

# -----------------------------------------------------------------------------
# SECTION 3 — Rg-derived PPFD as final fallback (note: applied in Section 6)
# -----------------------------------------------------------------------------
# This section documents the third fallback; the actual computation occurs in
# Section 6 after Rg is loaded from the Met Eirean solar file.
# Conversion: PPFD = Rg x 0.45 x 4.57
#   0.45 -- fraction of global radiation in the photosynthetically active range
#           (400-700 nm; McCree 1972)
#   4.57 -- conversion from W m-2 to umol photons m-2 s-1 in that waveband
# The formula is applied only where PPFD is still NA after Priorities 1 and 2,
# and only during daytime (Rg > 0); negative Rg at night is clipped to 0.

# -----------------------------------------------------------------------------
# SECTION 4 — Coalesce PPFD: Priority 1 (on-site) -> Priority 2 (twin-site)
# -----------------------------------------------------------------------------
# coalesce() fills NA values from PPFD_biomet with the corresponding twin-site
# value.  The final Priority 3 (Rg formula) coalesce is appended in Section 6
# after Rg is loaded.  Intermediate source columns are dropped to keep the
# data frame compact.
message("Coalescing on-site and twin-site PPFD ...")

twin_col <- paste0("PPFD_alt_", TWIN_SITE)

df <- df %>%
  mutate(PPFD = coalesce(PPFD_biomet, .data[[twin_col]])) %>%
  select(-PPFD_biomet, -any_of(twin_col))

message("  PPFD NA after on-site + twin coalesce: ", sum(is.na(df$PPFD)))

# -----------------------------------------------------------------------------
# SECTION 5 -- Met Eirean: Temp, RH, rain and derived VPD
# -----------------------------------------------------------------------------
# Met Eirean provides hourly observations; linear interpolation (zoo::na.approx,
# rule = 2) expands these to 30-min resolution.  rule = 2 extends the boundary
# values rather than producing NA at the edges of the year.
# VPD is derived from Temp and RH using the Tetens (Magnus) approximation:
#   e_s = (T_K^-8.2) * exp(77.345 + 0.0057*T_K - 7235/T_K)  [saturation VP, kPa]
#   e   = RH * e_s / 100                                      [actual VP, kPa]
#   VPD = e_s - e
# This formulation is consistent with the convention used in the REddyProc
# package and is numerically stable over the temperature range encountered.
message("Loading Met Eireann hourly weather ...")

met_path <- here::here("data", "meteireann_data", "met_hourly.csv")
if (!file.exists(met_path)) stop("Met Eireann file not found: ", met_path)

# Expected columns: date (format parseable by dmy_hm), temp, rhum, rain
met <- read_csv(met_path, show_col_types = FALSE) %>%
  mutate(timestamp = dmy_hm(date, tz = TZ)) %>%
  select(timestamp, Temp = temp, RH = rhum, rain) %>%
  filter(year(timestamp) == YEAR) %>%
  arrange(timestamp) %>%
  distinct(timestamp, .keep_all = TRUE)

df <- df %>%
  select(-RH) %>%
  left_join(met, by = "timestamp") %>%
  mutate(
    Temp = na.approx(Temp, x = timestamp, na.rm = FALSE, rule = 2),
    RH   = na.approx(RH,   x = timestamp, na.rm = FALSE, rule = 2),
    rain = na.approx(rain, x = timestamp, na.rm = FALSE, rule = 2),
    # VPD via Tetens (returns kPa)
    TempK = Temp + 273.15,
    e_s   = (TempK ^ -8.2) * exp(77.345 + 0.0057 * TempK - 7235 / TempK),
    e     = (RH * e_s) / 100,
    VPD   = e_s - e
  ) %>%
  select(-TempK, -e_s, -e)

# -----------------------------------------------------------------------------
# SECTION 6 -- Met Eirean: Global Radiation (Rg) and final PPFD fallback
# -----------------------------------------------------------------------------
# Hourly glorad (W m-2) from the solar file is linearly interpolated to 30 min.
# Rg = 0 is substituted for NA (missing radiation assumed zero, consistent with
# night-time or overcast periods where the measurement is near zero anyway).
# The Rg->PPFD formula (Section 3) is then applied and used as the Priority 3
# coalesce to fill any remaining PPFD gaps.
message("Loading Met Eireann solar radiation ...")

solar_path <- here::here("data", "meteireann_data", "solar_hourly.xlsx")
if (!file.exists(solar_path)) stop("Solar file not found: ", solar_path)

# Expected columns: date (POSIXct-parseable), glorad (W m⁻²)
solar <- read_excel(solar_path) %>%
  select(date, glorad) %>%
  mutate(
    timestamp = as.POSIXct(date, tz = TZ),
    glorad    = if_else(is.na(glorad), 0, as.numeric(glorad))
  ) %>%
  filter(year(timestamp) == YEAR) %>%
  arrange(timestamp) %>%
  distinct(timestamp, .keep_all = TRUE)

df <- df %>%
  left_join(select(solar, timestamp, glorad), by = "timestamp") %>%
  mutate(
    glorad = na.approx(glorad, x = timestamp, na.rm = FALSE, rule = 2),
    Rg     = glorad
  ) %>%
  select(-glorad)

# Rg → PPFD formula fills remaining gaps (PAR fraction 0.45, McCree factor 4.57)
# PPFD = Rg × 0.45 × 4.57  (applies only where PPFD is still NA after twin site)
df <- df %>%
  mutate(
    PPFD_from_Rg = pmax(Rg * 0.45 * 4.57, 0),   # negative Rg at night = 0
    PPFD         = coalesce(PPFD, PPFD_from_Rg)
  ) %>%
  select(-PPFD_from_Rg)

message("  PPFD NA after Rg formula fallback: ", sum(is.na(df$PPFD)))

# -----------------------------------------------------------------------------
# SECTION 7 -- Temporal encodings and season dummies
# -----------------------------------------------------------------------------
# Sine/cosine pairs encode cyclical time variables without discontinuities:
#   sin(2*pi*x/period) and cos(2*pi*x/period) jointly encode position in the
#   cycle so that, for example, 23:30 and 00:00 are close in predictor space.
# Season dummies use the Irish meteorological calendar:
#   Winter = Nov, Dec, Jan  |  Spring = Feb, Mar, Apr
#   Summer = May, Jun, Jul  |  Autumn = Aug, Sep, Oct
# This four-season encoding captures the broad seasonality of Irish grassland
# carbon exchange without imposing a single continuous seasonal curve.
message("Creating temporal predictors ...")

df <- df %>%
  mutate(
    .hr   = hour(timestamp) + minute(timestamp) / 60,
    .mon  = month(timestamp),
    .doy  = yday(timestamp),

    hour_sin  = sin(2 * pi * .hr  / 24),
    hour_cos  = cos(2 * pi * .hr  / 24),
    doy_sin   = sin(2 * pi * .doy / 365.25),
    doy_cos   = cos(2 * pi * .doy / 365.25),
    month_sin = sin(2 * pi * .mon / 12),
    month_cos = cos(2 * pi * .mon / 12),

    .season = case_when(
      .mon %in% c(12L,1L,2L) ~ "Winter",
      .mon %in% c(3L,4L,5L)  ~ "Spring",
      .mon %in% c(6L,7L,8L)  ~ "Summer",
      TRUE                   ~ "Autumn"
    ),

    Winter = as.integer(.season == "Winter"),
    Spring = as.integer(.season == "Spring"),
    Summer = as.integer(.season == "Summer"),
    Autumn = as.integer(.season == "Autumn")
  ) %>%
  select(-.hr, -.mon, -.doy, -.season)

# -----------------------------------------------------------------------------
# SECTION 8 -- Rolling rain (24-h trailing mean) and night flag
# -----------------------------------------------------------------------------
# rain_rolling_24 is a right-aligned (trailing) 48-half-hour rolling mean of
# rain.  partial = TRUE allows computation at the start of the year where fewer
# than 48 prior rows exist; fill = NA_real_ is used for the remaining edge cases.
# This predictor captures antecedent soil moisture conditions, which influence
# ecosystem respiration and stomatal conductance beyond the current half-hour.
#
# The night flag (night = 1 when PPFD < 10 umol m-2 s-1) is used internally
# by the flux-partitioning functions in the model scripts and by the night NEE
# filter applied here: night NEE < 0 is physiologically implausible and is
# set to NA to prevent the models from learning spurious negative night-time
# carbon uptake.
df <- df %>%
  arrange(timestamp) %>%
  mutate(
    rain_rolling_24 = zoo::rollapply(
      rain, width = 48,
      FUN   = function(z) mean(z, na.rm = TRUE),
      align = "right", partial = TRUE, fill = NA_real_
    ),
    night = if_else(is.finite(PPFD) & PPFD < 10, 1L, 0L, missing = 0L),
    # Night NEE < 0 is physiologically implausible
    NEE_orig = if_else(night == 1L & !is.na(NEE_orig) & NEE_orig < 0,
                        NA_real_, NEE_orig)
  )

# -----------------------------------------------------------------------------
# SECTION 9 -- Select final predictor columns and save
# -----------------------------------------------------------------------------
# Only the standardised predictor columns defined in core_preds are retained.
# Any extra columns added during intermediate computations (e.g. TempK, e_s)
# are dropped.  A warning is issued for any expected column that could not be
# constructed, allowing the issue to be diagnosed without stopping the pipeline.
core_preds <- c("timestamp","NEE_orig",
                "PPFD","Rg","VPD","RH","Temp","rain","rain_rolling_24","night",
                "hour_sin","hour_cos","doy_sin","doy_cos","month_sin","month_cos",
                "Winter","Spring","Summer","Autumn")

missing_cols <- setdiff(core_preds, names(df))
if (length(missing_cols))
  warning("Missing core columns: ", paste(missing_cols, collapse = ", "))

df_out <- df %>% select(all_of(intersect(core_preds, names(df))))

out_path <- file.path(OUT_DIR, paste0(SITE, "_", YEAR, "_predictors.rds"))
saveRDS(df_out, out_path)
message("Saved: ", out_path, "  (", nrow(df_out), " rows × ", ncol(df_out), " cols)")
message("\nScript 03 complete.")
