# =============================================================================
# 04_management_days_since.R -- Grazing_days_since & Fertiliser_days_since
# =============================================================================
#
# PURPOSE
#   Computes two continuous predictors measuring how long it has been since
#   the most recent management event of each type:
#
#     Grazing_days_since    -- days since the end of the most recent grazing event
#     Fertiliser_days_since -- days since the end of the most recent fertiliser event
#
# VALUE ENCODING
#   0        -- during an event window (the half-hour falls inside a grazing or
#               fertiliser event as defined by start_time/end_time in the CSV)
#   positive -- days elapsed since the end of the event
#
# RATIONALE FOR CONTINUOUS ENCODING
#   A continuous days-since variable captures the trajectory of ecosystem
#   recovery after disturbance.  Leaf area index, stomatal conductance, and
#   root respiration all change continuously during regrowth after grazing.
#   Encoding this as days-since allows the RF/MLP/XGBoost models to learn
#   nonlinear recovery curves from data, rather than imposing a binary
#   grazed/ungrazed classification that would discard recovery dynamics.
#
# YEAR-BOUNDARY CONTINUITY
#   The 'last event before year' lookup (last_events_before_year_{SITE}.csv)
#   provides the end timestamp of the most recent event that occurred before
#   the current year.  Without this, days_since at the start of January would
#   jump discontinuously when the two-year dataset is assembled in Script 07.
#
# FERTILISER DUMMY
#   A binary Fertiliser column (1 during fertiliser event windows, 0 otherwise)
#   is also added here as a by-product of the event-window processing.  It is
#   used by Script 05 to pin N step-function changes to the exact application
#   half-hour rather than midnight.
#
# INPUTS
#   data/data_prepared/{SITE}_{YEAR}_predictors.rds  (output of Script 03)
#
#   data/management_data/grazing_events_{SITE}.csv
#     Columns (one row per event day):
#       date        — YYYY-MM-DD
#       start_time  — HH:MM  (start of grazing window on that date)
#       end_time    — HH:MM or "24:00"
#       event_id    — integer grouping consecutive days in the same grazing bout
#
#   data/management_data/fertiliser_events_{SITE}.csv
#     Columns:
#       date        — YYYY-MM-DD
#       start_time  — HH:MM
#       end_time    — HH:MM
#       event_id    — integer
#
#   data/management_data/last_events_before_year_{SITE}.csv
#     Columns (one row per event type):
#       event_type  — "Grazing" | "Fertiliser"
#       last_end    — YYYY-MM-DD HH:MM:SS (UTC)  of the last event end before
#                     the start of the current YEAR
#     Purpose: back-fill days_since at the start of the year when no event
#     occurred yet in the current YEAR.
#
# OUTPUTS
#   data/data_prepared/{SITE}_{YEAR}_predictors.rds  (updated in-place with
#     Grazing_days_since and Fertiliser_days_since added)
# =============================================================================

library(here)
library(dplyr)
library(readr)
library(lubridate)
library(stringr)
library(purrr)

# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------
TZ   <- "UTC"

# -----------------------------------------------------------------------------
# SECTION 1 — Load predictors RDS
# -----------------------------------------------------------------------------
in_path <- here::here("data", "data_prepared",
                       paste0(SITE, "_", YEAR, "_predictors.rds"))
if (!file.exists(in_path)) stop("Predictors RDS not found: ", in_path)
df <- readRDS(in_path)
df <- df %>%
  mutate(timestamp = as.POSIXct(timestamp, tz = TZ)) %>%
  arrange(timestamp)
message("Loaded: ", in_path, "  (", nrow(df), " rows)")

# -----------------------------------------------------------------------------
# HELPER: compute_days_since()
# -----------------------------------------------------------------------------
# Core function called twice: once for grazing events, once for fertiliser.
# Algorithm:
#   1. If no events exist for this year, return days since last_end_prev for
#      all rows (monotonically increasing from the previous year's last event).
#   2. Parse start and end times from the event CSV.  end_time == '24:00' is
#      converted to 00:00 of the following day (representing end of day).
#   3. Collapse all rows sharing the same event_id into a single time window
#      (min start -> max end) to handle multi-day events stored as daily rows.
#   4. For each half-hour timestamp:
#      a. If inside any event window -> days_since = 0.
#      b. If after the most recent event end -> days_since = difftime(ts, max_end).
#      c. If before any event this year -> days_since = difftime(ts, last_end_prev).
#   5. Final pmax(days_vec, 0) prevents negative values from timing edge cases.
# Arguments:
#   df_ts         -- data frame with column 'timestamp'
#   events        -- data frame: date, start_time, end_time, event_id
#   last_end_prev -- POSIXct: last event end before the current year
# Returns: numeric vector, length = nrow(df_ts)
# -----------------------------------------------------------------------------
compute_days_since <- function(df_ts,            # tibble with column 'timestamp'
                                events,           # tibble: date, start_time, end_time
                                last_end_prev,    # POSIXct: last event end before year
                                tz = "UTC") {
  n <- nrow(df_ts)
  days_vec <- rep(NA_real_, n)

  if (nrow(events) == 0) {
    # No events in this year: entire year is days-since previous event
    days_vec <- as.numeric(difftime(df_ts$timestamp, last_end_prev, units = "days"))
    days_vec <- pmax(days_vec, 0)
    return(days_vec)
  }

  # Build POSIX windows for every event row
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
    # Collapse all windows of same event_id into one span (min start → max end)
    group_by(event_id) %>%
    summarise(start_ct = min(start_ct), end_ct = max(end_ct), .groups = "drop") %>%
    arrange(start_ct)

  # For each timestamp, compute days since nearest preceding event end
  # or set to 0 if inside an event window
  event_ends <- sort(windows$end_ct)

  for (i in seq_len(n)) {
    ts <- df_ts$timestamp[i]

    # Is ts inside ANY event window?
    inside <- any(ts >= windows$start_ct & ts < windows$end_ct)
    if (inside) {
      days_vec[i] <- 0
      next
    }

    # Find the most recent event end before ts
    past_ends <- event_ends[event_ends <= ts]
    if (length(past_ends) > 0) {
      days_vec[i] <- as.numeric(difftime(ts, max(past_ends), units = "days"))
    } else {
      # No event yet in this year: use previous year's last event
      if (!is.null(last_end_prev) && !is.na(last_end_prev)) {
        days_vec[i] <- as.numeric(difftime(ts, last_end_prev, units = "days"))
      }
    }
  }

  pmax(days_vec, 0, na.rm = FALSE)
}

# -----------------------------------------------------------------------------
# SECTION 2 -- Load event tables
# -----------------------------------------------------------------------------
# load_events() reads the per-site grazing or fertiliser event CSV and filters
# to the current year.  Returns an empty tibble if the file is absent, allowing
# the script to continue without management data for years with no events.
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

# Load previous-year last event ends
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

message("Grazing events in ", YEAR, ": ", nrow(graz_events), " rows")
message("Fertiliser events in ", YEAR, ": ", nrow(fert_events), " rows")
message("Last grazing end before ", YEAR, ": ", last_graz_prev)
message("Last fertiliser end before ", YEAR, ": ", last_fert_prev)

# -----------------------------------------------------------------------------
# SECTION 3 -- Compute days_since and Fertiliser dummy
# -----------------------------------------------------------------------------
# Calls compute_days_since() for both event types.
# Also constructs the Fertiliser binary dummy column: 1 = the current half-hour
# falls inside a fertiliser event window, 0 = outside.  Event windows are
# derived from the same fert_events table used for days_since, without reference
# to the event_id (each row is its own half-hourly window check).
message("Computing Grazing_days_since ...")
df$Grazing_days_since <- compute_days_since(
  df, graz_events, last_graz_prev, tz = TZ
)

message("Computing Fertiliser_days_since ...")
df$Fertiliser_days_since <- compute_days_since(
  df, fert_events, last_fert_prev, tz = TZ
)

message("  Grazing_days_since range: [",
        round(min(df$Grazing_days_since, na.rm=TRUE), 1), ", ",
        round(max(df$Grazing_days_since, na.rm=TRUE), 1), "]")
message("  Fertiliser_days_since range: [",
        round(min(df$Fertiliser_days_since, na.rm=TRUE), 1), ", ",
        round(max(df$Fertiliser_days_since, na.rm=TRUE), 1), "]")

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
message("  Fertiliser dummy: ", sum(df$Fertiliser), " half-hours inside event windows")
# -----------------------------------------------------------------------------
# SECTION 4 -- Save updated RDS (in-place)
# -----------------------------------------------------------------------------
# Overwrites the same RDS file produced by Script 03, adding three new columns:
#   Grazing_days_since, Fertiliser_days_since, Fertiliser (dummy).
# All other columns are unchanged.
saveRDS(df, in_path)
message("Updated RDS saved: ", in_path)
message("\nScript 04 complete.")
