# 03_predictors_meteo.R — meteorological predictors, gap-filling & time encodings
#
# Builds the BASE predictor set for one site x year from the full-grid QC data.
# Two branches share all downstream logic (temporal encodings, VPD, rolling
# rain, night flag):
#   JC1 2020 (IS_OPEN_PATH == TRUE):
#     NEE_orig = co2_flux_base_filters_3_70_grid (open-path column)
#     PPFD / Temp taken from on-site biomet as-is; RH, rain, Rg from Met Éireann
#   2023 / 2024 (IS_OPEN_PATH == FALSE):
#     NEE_orig = co2_flux_base_filters_6_70_grid (closed-path column)
#     PPFD gap-fill cascade on-site -> twin-site -> Rg formula; Met Éireann meteo
#
# Input  : data/data_qc/results_qc_{SITE}_{YEAR}/merged_qc_fullgrid.csv
# Output : data/data_prepared/{SITE}_{YEAR}_predictors.rds

library(here)
library(dplyr)
library(readr)
library(readxl)
library(lubridate)
library(tidyr)
library(zoo)
library(ranger)

# settings (SITE, YEAR, IS_OPEN_PATH, TWIN_SITE injected by script 01)
if (!exists("SITE"))         SITE         <- "JC1"
if (!exists("YEAR"))         YEAR         <- 2023
if (!exists("IS_OPEN_PATH")) IS_OPEN_PATH <- FALSE
if (!exists("TWIN_SITE"))    TWIN_SITE    <- setdiff(c("JC1", "JC2"), SITE)
TZ <- "UTC"

PPFD_BADPERIOD_DAYS  <- 3
PPFD_ALLLT_THRESHOLD <- 375
PPFD_ALLLT_DAYS      <- 7

OUT_DIR <- here::here("data", "data_prepared")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# flag contiguous runs of implausible PPFD (out-of-range, or persistently low)
detect_bad_periods <- function(df, col, days = 3L, interval_per_day = 48L,
                               all_lt_threshold = 375, all_lt_days = 7L) {
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
  if (!nrow(segs))
    return(tibble::tibble(start_idx = integer(), end_idx = integer(),
                          start_time = as.POSIXct(character()),
                          end_time   = as.POSIXct(character())))
  
  merged_segs <- list(); cur <- segs[1L, ]
  for (i in seq_len(nrow(segs))[-1]) {
    nxt <- segs[i, ]
    if (nxt$start_idx <= cur$end_idx + 1L)
      cur$end_idx <- max(cur$end_idx, nxt$end_idx)
    else { merged_segs[[length(merged_segs) + 1L]] <- cur; cur <- nxt }
  }
  merged_segs[[length(merged_segs) + 1L]] <- cur
  dplyr::bind_rows(merged_segs) %>%
    dplyr::mutate(start_time = df$timestamp[start_idx],
                  end_time   = df$timestamp[end_idx])
}

# set the flagged bad-period rows of a column to NA
apply_bad_periods_na <- function(df, col, segs) {
  if (!nrow(segs)) return(df)
  bad <- rep(FALSE, nrow(df))
  for (i in seq_len(nrow(segs)))
    bad[segs$start_idx[i]:segs$end_idx[i]] <- TRUE
  df[[col]][bad] <- NA_real_
  df
}

# load the full-grid QC data for this site-year
in_path <- here::here("data", "data_qc",
                      paste0("results_qc_", SITE, "_", YEAR),
                      "merged_qc_fullgrid.csv")
if (!file.exists(in_path)) stop("Input not found: ", in_path)

df <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(timestamp = as.POSIXct(timestamp, tz = TZ)) %>%
  arrange(timestamp) %>%
  distinct(timestamp, .keep_all = TRUE)

if (IS_OPEN_PATH) {
  # JC1 2020 open-path branch
  
  # NEE: open-path QC column (stationarity filter 3, footprint >= 70%)
  if (!"NEE_orig" %in% names(df)) {
    if ("co2_flux_base_filters_3_70_grid" %in% names(df)) {
      df <- df %>% rename(NEE_orig = co2_flux_base_filters_3_70_grid)
    } else {
      stop("co2_flux_base_filters_3_70_grid not found in 2020 QC grid file.")
    }
  }
  
  # PPFD: on-site biomet used as-is (bad-period filter still applied)
  if ("PPFD_biomet" %in% names(df)) {
    df <- df %>%
      mutate(PPFD_biomet = na_if(PPFD_biomet, -9999),
             PPFD_biomet = if_else(!is.na(PPFD_biomet) & PPFD_biomet < 0,
                                   -PPFD_biomet, PPFD_biomet))
    bp <- detect_bad_periods(df, "PPFD_biomet",
                             days           = PPFD_BADPERIOD_DAYS,
                             all_lt_threshold = PPFD_ALLLT_THRESHOLD,
                             all_lt_days    = PPFD_ALLLT_DAYS)
    df <- apply_bad_periods_na(df, "PPFD_biomet", bp)
    df <- df %>% rename(PPFD = PPFD_biomet)
  } else {
    warning("PPFD_biomet not found in 2020 QC file — PPFD set to NA.")
    df$PPFD <- NA_real_
  }
  
  # Temp: on-site biomet air_temperature (convert from Kelvin if needed)
  if ("air_temperature" %in% names(df)) {
    df <- df %>%
      mutate(Temp = if_else(!is.na(air_temperature) & air_temperature > 200,
                            air_temperature - 273.15,
                            air_temperature))
  } else {
    warning("air_temperature not in 2020 QC file — Temp set to NA.")
    df$Temp <- NA_real_
  }
  
  # RH & rain: Met Éireann hourly (on-site RH may have gaps in 2020)
  met_path <- here::here("data", "meteireann_data", "met_hourly.csv")
  if (!file.exists(met_path)) stop("Met Eireann file not found: ", met_path)
  
  met <- read_csv(met_path, show_col_types = FALSE) %>%
    mutate(timestamp = dmy_hm(date, tz = TZ)) %>%
    select(timestamp, RH = rhum, rain) %>%
    filter(year(timestamp) == YEAR) %>%
    arrange(timestamp) %>%
    distinct(timestamp, .keep_all = TRUE)
  
  # drop on-site RH (replaced by Met Éireann), then join and interpolate
  df <- df %>% select(-any_of("RH"))
  
  df <- df %>%
    left_join(met, by = "timestamp") %>%
    mutate(
      RH   = na.approx(RH,   x = timestamp, na.rm = FALSE, rule = 2),
      rain = na.approx(rain, x = timestamp, na.rm = FALSE, rule = 2)
    )
  
  # VPD from Temp + RH
  df <- df %>%
    mutate(
      TempK = Temp + 273.15,
      e_s   = (TempK ^ -8.2) * exp(77.345 + 0.0057 * TempK - 7235 / TempK),
      e     = (RH * e_s) / 100,
      VPD   = e_s - e
    ) %>%
    select(-TempK, -e_s, -e)
  
  # Rg: Rg_biomet already in the QC file — gap-fill the few missing points
  df <- df %>%
    mutate(
      Rg = na_if(as.numeric(Rg_biomet), -9999),
      Rg = na.approx(Rg, x = timestamp, na.rm = FALSE, rule = 2),
      Rg = pmax(Rg, 0, na.rm = FALSE)
    )
  
} else {
  # standard closed-path branch (2023 / 2024)
  
  # NEE: closed-path QC column (stationarity filter 6, footprint >= 70%)
  if (!"NEE_orig" %in% names(df) && "co2_flux_base_filters_6_70_grid" %in% names(df))
    df <- df %>% rename(NEE_orig = co2_flux_base_filters_6_70_grid)
  
  # PPFD: on-site biomet, sign-corrected and bad-period filtered
  ppfd_biomet_col <- intersect(c("PPFD_1_1_1_biomet", "PPFD_1_1_1", "PPFD_biomet"),
                               names(df))[1]
  if (is.na(ppfd_biomet_col)) {
    df$PPFD_biomet <- NA_real_
  } else {
    df <- df %>%
      rename(PPFD_biomet = all_of(ppfd_biomet_col)) %>%
      mutate(PPFD_biomet = na_if(PPFD_biomet, -9999),
             PPFD_biomet = if_else(!is.na(PPFD_biomet) & PPFD_biomet < 0,
                                   -PPFD_biomet, PPFD_biomet))
  }
  
  bp <- detect_bad_periods(df, "PPFD_biomet",
                           days           = PPFD_BADPERIOD_DAYS,
                           all_lt_threshold = PPFD_ALLLT_THRESHOLD,
                           all_lt_days    = PPFD_ALLLT_DAYS)
  df <- apply_bad_periods_na(df, "PPFD_biomet", bp)
  
  # twin-site PPFD fallback (the other grassland site's biomet)
  twin_ppfd_col    <- paste0("PPFD_alt_", TWIN_SITE)
  twin_biomet_path <- here::here("data", "raw_data",
                                 paste0(TWIN_SITE, "_", YEAR, "_biomet.xlsx"))
  
  if (file.exists(twin_biomet_path)) {
    twin_df <- read_excel(twin_biomet_path) %>%
      distinct(date, time, .keep_all = TRUE) %>%
      mutate(timestamp = update(ymd_hms(time, tz = "Europe/Dublin"),
                                year  = year(as_date(date)),
                                month = month(as_date(date)),
                                mday  = mday(as_date(date)))) %>%
      arrange(timestamp) %>%
      distinct(timestamp, .keep_all = TRUE) %>%
      select(timestamp, any_of(c("PPFD_1_1_1_biomet", "PPFD_1_1_1"))) %>%
      rename_with(~ twin_ppfd_col, -timestamp) %>%
      mutate(!!twin_ppfd_col := na_if(.data[[twin_ppfd_col]], -9999),
             !!twin_ppfd_col := if_else(!is.na(.data[[twin_ppfd_col]]) &
                                          .data[[twin_ppfd_col]] < 0,
                                        -.data[[twin_ppfd_col]],
                                        .data[[twin_ppfd_col]]))
    df <- df %>% left_join(twin_df, by = "timestamp")
    bp2 <- detect_bad_periods(df, twin_ppfd_col,
                              days           = PPFD_BADPERIOD_DAYS,
                              all_lt_threshold = PPFD_ALLLT_THRESHOLD,
                              all_lt_days    = PPFD_ALLLT_DAYS)
    df <- apply_bad_periods_na(df, twin_ppfd_col, bp2)
  } else {
    df[[twin_ppfd_col]] <- NA_real_
  }
  
  # combine on-site and twin PPFD (on-site takes priority)
  df <- df %>%
    mutate(PPFD = coalesce(PPFD_biomet, .data[[twin_ppfd_col]])) %>%
    select(-PPFD_biomet, -any_of(twin_ppfd_col))
  
  # Met Éireann Temp, RH, rain, then VPD from Temp + RH
  met_path <- here::here("data", "meteireann_data", "met_hourly.csv")
  if (!file.exists(met_path)) stop("Met Eireann file not found: ", met_path)
  
  met <- read_csv(met_path, show_col_types = FALSE) %>%
    mutate(timestamp = dmy_hm(date, tz = TZ)) %>%
    select(timestamp, Temp = temp, RH = rhum, rain) %>%
    filter(year(timestamp) == YEAR) %>%
    arrange(timestamp) %>%
    distinct(timestamp, .keep_all = TRUE)
  
  df <- df %>%
    select(-any_of("RH")) %>%
    left_join(met, by = "timestamp") %>%
    mutate(
      Temp = na.approx(Temp, x = timestamp, na.rm = FALSE, rule = 2),
      RH   = na.approx(RH,   x = timestamp, na.rm = FALSE, rule = 2),
      rain = na.approx(rain, x = timestamp, na.rm = FALSE, rule = 2),
      TempK = Temp + 273.15,
      e_s   = (TempK ^ -8.2) * exp(77.345 + 0.0057 * TempK - 7235 / TempK),
      e     = (RH * e_s) / 100,
      VPD   = e_s - e
    ) %>%
    select(-TempK, -e_s, -e)
  
  # Met Éireann global radiation -> Rg
  solar_path <- here::here("data", "meteireann_data", "solar_hourly.xlsx")
  if (!file.exists(solar_path)) stop("Solar file not found: ", solar_path)
  
  solar <- read_excel(solar_path) %>%
    select(date, glorad) %>%
    mutate(timestamp = as.POSIXct(date, tz = TZ),
           glorad    = if_else(is.na(glorad), 0, as.numeric(glorad))) %>%
    filter(year(timestamp) == YEAR) %>%
    arrange(timestamp) %>%
    distinct(timestamp, .keep_all = TRUE)
  
  df <- df %>%
    left_join(select(solar, timestamp, glorad), by = "timestamp") %>%
    mutate(glorad = na.approx(glorad, x = timestamp, na.rm = FALSE, rule = 2),
           Rg     = glorad) %>%
    select(-glorad)
  
  # fill any remaining PPFD gaps from Rg (PPFD = Rg x 0.45 x 4.57)
  df <- df %>%
    mutate(PPFD_from_Rg = pmax(Rg * 0.45 * 4.57, 0),
           PPFD         = coalesce(PPFD, PPFD_from_Rg)) %>%
    select(-PPFD_from_Rg)
}

# temporal predictors: cyclical hour/doy/month encodings + season dummies
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
      .mon %in% c(12L, 1L, 2L) ~ "Winter",
      .mon %in% c(3L,  4L, 5L) ~ "Spring",
      .mon %in% c(6L,  7L, 8L) ~ "Summer",
      TRUE                     ~ "Autumn"
    ),
    Winter = as.integer(.season == "Winter"),
    Spring = as.integer(.season == "Spring"),
    Summer = as.integer(.season == "Summer"),
    Autumn = as.integer(.season == "Autumn")
  ) %>%
  select(-.hr, -.mon, -.doy, -.season)

# 24 h rolling rain, night flag, and drop negative night-time NEE
df <- df %>%
  arrange(timestamp) %>%
  mutate(
    rain_rolling_24 = zoo::rollapply(
      rain, width = 48,
      FUN   = function(z) mean(z, na.rm = TRUE),
      align = "right", partial = TRUE, fill = NA_real_
    ),
    night    = if_else(is.finite(PPFD) & PPFD < 10, 1L, 0L, missing = 0L),
    NEE_orig = if_else(night == 1L & !is.na(NEE_orig) & NEE_orig < 0,
                       NA_real_, NEE_orig)
  )

# select the BASE predictor columns and save
core_preds <- c("timestamp", "NEE_orig",
                "PPFD", "Rg", "VPD", "RH", "Temp", "rain", "rain_rolling_24",
                "night",
                "hour_sin", "hour_cos", "doy_sin", "doy_cos",
                "month_sin", "month_cos",
                "Winter", "Spring", "Summer", "Autumn")

missing_cols <- setdiff(core_preds, names(df))
if (length(missing_cols))
  warning("Missing core columns: ", paste(missing_cols, collapse = ", "))

df_out <- df %>% select(all_of(intersect(core_preds, names(df))))

out_path <- file.path(OUT_DIR, paste0(SITE, "_", YEAR, "_predictors.rds"))
saveRDS(df_out, out_path)
