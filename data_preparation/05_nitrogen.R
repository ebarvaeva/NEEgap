# =============================================================================
# 05_nitrogen.R -- Nitrogen Application Step Function
# =============================================================================
#
# PURPOSE
#   Creates the predictor N (kg N ha-1): the nitrogen amount applied at the
#   most recent fertiliser event.  The value is a step function:
#
#     N = 0            before the first fertiliser event of the year
#     N = N_kg_ha[i]   at and after the i-th event, until the next event
#
# WHY A STEP FUNCTION
#   Unlike temperature or radiation, N availability in the soil does not
#   reset instantly between applications.  A step function encoding captures
#   the 'current N load' at the time of each gap, allowing the model to learn
#   whether a given NEE observation occurred under high or low N fertilisation.
#   The step persists until superseded by the next application, reflecting
#   the slow mineralisation of urea and calcium ammonium nitrate in the soil.
#
# SLURRY EXCLUSION
#   Only mineral-N events (is_slurry == FALSE) are included.  Slurry
#   applications are excluded because their effective N content is uncertain
#   (it depends on dry matter content and application method) and is not
#   directly comparable to the certified N content of mineral fertilisers.
#
# APPLICATION TIME ALIGNMENT
#   The N step uses the first half-hour where the Fertiliser dummy (from
#   Script 04) equals 1 on the event date, rather than midnight.  This aligns
#   the step with the actual application time as recorded in the event CSV.
#   If the Fertiliser dummy is absent, midnight of the event date is used.
#
# INPUTS
#   data/data_prepared/{SITE}_{YEAR}_predictors.rds  (must already have
#     Fertiliser_days_since from Script 04)
#
#   data/management_data/nitrogen_amounts_{SITE}.csv
#     Columns:
#       fert_date   — DD/MM/YYYY  date of fertiliser application
#       N_kg_ha     — numeric, kg N ha⁻¹ applied on that date
#       is_slurry   — logical (TRUE/FALSE); slurry rows are excluded
#
# OUTPUTS
#   data/data_prepared/{SITE}_{YEAR}_predictors.rds  (updated in-place with N)
#
# LOGIC
#   N is assigned as a step function using findInterval():
#     - Before the first event of the year: N = 0
#     - At/after event i: N = N_kg_ha[i]  (value stays until next event)
#   The "Fertiliser == 1" dummy (built in Script 04 from the event windows) is
#   used to locate the exact first observed half-hour of each event, so the
#   step aligns with the actual application time rather than midnight.
# =============================================================================

library(here)
library(dplyr)
library(readr)
library(lubridate)
library(stringr)

# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------
TZ   <- "UTC"

# -----------------------------------------------------------------------------
# SECTION 1 -- Load predictors RDS
# -----------------------------------------------------------------------------
# Loads the predictors RDS produced and updated by Scripts 03-04.
# The Fertiliser dummy column (added by Script 04) is checked here;
# if absent, a warning is issued and N step changes default to midnight.
in_path <- here::here("data", "data_prepared",
                       paste0(SITE, "_", YEAR, "_predictors.rds"))
if (!file.exists(in_path)) stop("Predictors RDS not found: ", in_path)
df <- readRDS(in_path) %>%
  mutate(timestamp = as.POSIXct(timestamp, tz = TZ)) %>%
  arrange(timestamp)
message("Loaded: ", in_path)

# Fertiliser dummy is a by-product of Script 04's event-window processing.
# If absent (e.g. Script 04 was not run), application time defaults to midnight.
has_fert_dummy <- "Fertiliser" %in% names(df)
if (!has_fert_dummy)
  message("  'Fertiliser' dummy not found — application time will default to midnight.")

# -----------------------------------------------------------------------------
# SECTION 2 -- Load nitrogen amounts
# -----------------------------------------------------------------------------
# Reads nitrogen_amounts_{SITE}.csv, excludes slurry rows, and filters to
# events up to the last observed timestamp in df.  Events from future years
# are excluded because the step function must not anticipate future applications.
n_path <- here::here("data", "management_data",
                      paste0("nitrogen_amounts_", SITE, ".csv"))
if (!file.exists(n_path)) stop("Nitrogen amounts file not found: ", n_path)

n_tbl <- read_csv(n_path, show_col_types = FALSE) %>%
  mutate(fert_date = dmy(str_replace_all(fert_date, "\\.", "/"))) %>%
  filter(if ("is_slurry" %in% names(.)) !is_slurry else TRUE) %>%
  filter(fert_date <= max(df$timestamp)) %>%
  arrange(fert_date)

message("Nitrogen events in ", YEAR, ": ", nrow(n_tbl))
if (nrow(n_tbl)) print(n_tbl %>% select(fert_date, N_kg_ha))

# -----------------------------------------------------------------------------
# SECTION 3 -- Compute N step function
# -----------------------------------------------------------------------------
# findInterval(df$timestamp, n_starts) returns index i such that
#   n_starts[i] <= timestamp < n_starts[i+1]
# Index 0 means the timestamp precedes all events -> N = 0.
# Index i > 0 means the most recent event was the i-th one -> N = n_vals[i].
# pmax(idx, 1L) prevents out-of-bounds indexing; the idx == 0 case is handled
# by the if_else before indexing into n_vals.
if (nrow(n_tbl) == 0) {
  # No nitrogen events this year
  df$N <- 0
} else {
  # Locate the first Fertiliser == 1 half-hour on each event date, or default
  # to midnight if the dummy is missing
  if (has_fert_dummy) {
    fert_starts <- df %>%
      filter(Fertiliser == 1) %>%
      mutate(fert_date = as.Date(timestamp)) %>%
      group_by(fert_date) %>%
      summarise(start_dt = min(timestamp), .groups = "drop")
  } else {
    fert_starts <- tibble(
      fert_date = as.Date(character()),
      start_dt  = as.POSIXct(character(), tz = TZ)
    )
  }

  fert_schedule <- n_tbl %>%
    left_join(fert_starts, by = "fert_date") %>%
    mutate(
      start_dt = coalesce(start_dt, as.POSIXct(fert_date, tz = TZ))
    ) %>%
    arrange(start_dt)

  n_starts <- fert_schedule$start_dt
  n_vals   <- fert_schedule$N_kg_ha

  # findInterval returns index i such that n_starts[i] <= timestamp < n_starts[i+1]
  # index 0 means before the first event → N = 0
  idx <- findInterval(df$timestamp, n_starts)
  df$N <- if_else(idx == 0L, 0, as.numeric(n_vals[pmax(idx, 1L)]))
}

message("  N range: [", min(df$N, na.rm=TRUE), ", ", max(df$N, na.rm=TRUE), "]")

# -----------------------------------------------------------------------------
# SECTION 4 -- Save updated RDS
# -----------------------------------------------------------------------------
# Overwrites the per-year predictors RDS, adding column N.
saveRDS(df, in_path)
message("Updated RDS saved: ", in_path)
message("\nScript 05 complete.")
