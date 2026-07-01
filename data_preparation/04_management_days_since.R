# 04_management_days_since.R — Grazing_days_since & Fertiliser_days_since
#
# Adds two continuous predictors measuring days elapsed since the most recent
# grazing / fertiliser event ended: 0 during an event window, positive after it.
# Start-of-year values are back-filled from the previous year's last event so
# 2023 and 2024 join seamlessly when script 07 concatenates them. Also adds the
# Fertiliser 0/1 dummy used by script 05 to pin the N step to the event time.
#
# Inputs :
#   data/data_prepared/{SITE}_{YEAR}_predictors.rds            (from script 03)
#   data/management_data/grazing_events_{SITE}.csv             (date, start_time, end_time, event_id)
#   data/management_data/fertiliser_events_{SITE}.csv          (same columns)
#   data/management_data/last_events_before_year_{SITE}.csv    (event_type, last_end per year)
# Output :
#   data/data_prepared/{SITE}_{YEAR}_predictors.rds  (updated in-place)

library(here)
library(dplyr)
library(readr)
library(lubridate)
library(stringr)
library(purrr)

TZ   <- "UTC"

# load the predictors RDS for this site-year
in_path <- here::here("data", "data_prepared",
                      paste0(SITE, "_", YEAR, "_predictors.rds"))
if (!file.exists(in_path)) stop("Predictors RDS not found: ", in_path)
df <- readRDS(in_path)
df <- df %>%
  mutate(timestamp = as.POSIXct(timestamp, tz = TZ)) %>%
  arrange(timestamp)

# days_since = 0 inside event windows, positive days after end, back-filled
# from last_end_prev before the first event of the current year
compute_days_since <- function(df_ts,            # tibble with column 'timestamp'
                               events,           # tibble: date, start_time, end_time
                               last_end_prev,    # POSIXct: last event end before year
                               tz = "UTC") {
  n <- nrow(df_ts)
  days_vec <- rep(NA_real_, n)
  
  if (nrow(events) == 0) {
    # no events this year: whole year is days-since the previous event
    days_vec <- as.numeric(difftime(df_ts$timestamp, last_end_prev, units = "days"))
    days_vec <- pmax(days_vec, 0)
    return(days_vec)
  }
  
  # collapse each event_id into one span (min start -> max end)
  windows <- events %>%
    mutate(
      date  = as.Date(date),
      end_is_2400 = end_time == "24:00",
      start_ct = as.POSIXct(paste(date, str_pad(start_time, 5, pad="0")), tz = tz),
      end_ct   = if_else(
        end_is_2400,
        as.POSIXct(paste(date + days(1), "00:00"), tz = tz),
        as.POSIXct(paste(date, str_pad(end_time, 5, pad="0")), tz = tz)
      )
    ) %>%
    filter(!is.na(start_ct), !is.na(end_ct), end_ct > start_ct) %>%
    group_by(event_id) %>%
    summarise(start_ct = min(start_ct), end_ct = max(end_ct), .groups = "drop") %>%
    arrange(start_ct)
  
  event_ends <- sort(windows$end_ct)
  
  # for each timestamp: 0 if inside a window, else days since nearest past end
  for (i in seq_len(n)) {
    ts <- df_ts$timestamp[i]
    
    inside <- any(ts >= windows$start_ct & ts < windows$end_ct)
    if (inside) {
      days_vec[i] <- 0
      next
    }
    
    past_ends <- event_ends[event_ends <= ts]
    if (length(past_ends) > 0) {
      days_vec[i] <- as.numeric(difftime(ts, max(past_ends), units = "days"))
    } else {
      # no event yet this year: use previous year's last event end
      if (!is.null(last_end_prev) && !is.na(last_end_prev)) {
        days_vec[i] <- as.numeric(difftime(ts, last_end_prev, units = "days"))
      }
    }
  }
  
  pmax(days_vec, 0, na.rm = FALSE)
}

# read one event-type CSV, filtered to the current year
load_events <- function(path, year, tz) {
  if (!file.exists(path)) {
    warning("Events file not found: ", path); return(tibble())
  }
  read_csv(path, show_col_types = FALSE) %>%
    mutate(date = as.Date(date)) %>%
    filter(year(date) == year) %>%
    arrange(date, start_time)
}

path_graz <- here::here("data", "management_data",
                        paste0("grazing_events_", SITE, ".csv"))
path_fert <- here::here("data", "management_data",
                        paste0("fertiliser_events_", SITE, ".csv"))
path_prev <- here::here("data", "management_data",
                        paste0("last_events_before_year_", SITE, ".csv"))

graz_events <- load_events(path_graz, YEAR, TZ)
fert_events <- load_events(path_fert, YEAR, TZ)

# previous-year last event ends (fallback: Jan 1 of the current year)
if (file.exists(path_prev)) {
  prev_ends <- read_csv(path_prev, show_col_types = FALSE) %>%
    filter(year_for == YEAR) %>%
    mutate(last_end = as.POSIXct(last_end, tz = TZ))
  last_graz_prev <- prev_ends %>%
    filter(event_type == "Grazing") %>% pull(last_end) %>% first()
  last_fert_prev <- prev_ends %>%
    filter(event_type == "Fertiliser") %>% pull(last_end) %>% first()
} else {
  warning("last_events_before_year file not found — using Jan 1 of year as fallback.")
  last_graz_prev <- as.POSIXct(paste0(YEAR, "-01-01"), tz = TZ)
  last_fert_prev <- as.POSIXct(paste0(YEAR, "-01-01"), tz = TZ)
}

# days-since columns for grazing and fertiliser
df$Grazing_days_since <- compute_days_since(
  df, graz_events, last_graz_prev, tz = TZ
)

df$Fertiliser_days_since <- compute_days_since(
  df, fert_events, last_fert_prev, tz = TZ
)

# Fertiliser 0/1 dummy: 1 inside each fertiliser window (used by script 05)
fert_windows <- fert_events %>%
  mutate(
    end_is_2400 = end_time == "24:00",
    start_ct = as.POSIXct(paste(as.Date(date), str_pad(start_time, 5, pad="0")), tz = TZ),
    end_ct   = if_else(
      end_is_2400,
      as.POSIXct(paste(as.Date(date) + days(1), "00:00"), tz = TZ),
      as.POSIXct(paste(as.Date(date), str_pad(end_time, 5, pad="0")), tz = TZ)
    )
  ) %>%
  filter(!is.na(start_ct), end_ct > start_ct)

df$Fertiliser <- 0L
for (i in seq_len(nrow(fert_windows))) {
  inside <- df$timestamp >= fert_windows$start_ct[i] &
    df$timestamp <  fert_windows$end_ct[i]
  df$Fertiliser[inside] <- 1L
}

# save updated RDS in place
saveRDS(df, in_path)
