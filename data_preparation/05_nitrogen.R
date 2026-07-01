# 05_nitrogen.R — nitrogen applied (N, kg N ha-1)
#
# Adds N as a step function: it jumps to the amount applied at each mineral
# fertiliser event and holds until the next event (0 before the first event).
# Slurry rows are excluded (is_slurry == TRUE). The Fertiliser dummy from
# script 04 pins each step to the first observed half-hour of the event
# (falls back to midnight if the dummy is absent). Assigned with findInterval().
#
# Inputs :
#   data/data_prepared/{SITE}_{YEAR}_predictors.rds        (must have Fertiliser_days_since)
#   data/management_data/nitrogen_amounts_{SITE}.csv       (fert_date, N_kg_ha, is_slurry)
# Output :
#   data/data_prepared/{SITE}_{YEAR}_predictors.rds  (updated in-place with N)

library(here)
library(dplyr)
library(readr)
library(lubridate)
library(stringr)

TZ   <- "UTC"

# load the predictors RDS for this site-year
in_path <- here::here("data", "data_prepared",
                      paste0(SITE, "_", YEAR, "_predictors.rds"))
if (!file.exists(in_path)) stop("Predictors RDS not found: ", in_path)
df <- readRDS(in_path) %>%
  mutate(timestamp = as.POSIXct(timestamp, tz = TZ)) %>%
  arrange(timestamp)

# the Fertiliser dummy (script 04) pins the step to the event's first half-hour
has_fert_dummy <- "Fertiliser" %in% names(df)

# load nitrogen amounts, dropping slurry rows and events past the data window
n_path <- here::here("data", "management_data",
                     paste0("nitrogen_amounts_", SITE, ".csv"))
if (!file.exists(n_path)) stop("Nitrogen amounts file not found: ", n_path)

n_tbl <- read_csv(n_path, show_col_types = FALSE) %>%
  mutate(fert_date = dmy(str_replace_all(fert_date, "\\.", "/"))) %>%
  filter(if ("is_slurry" %in% names(.)) !is_slurry else TRUE) %>%
  filter(fert_date <= max(df$timestamp)) %>%
  arrange(fert_date)

# N step function via findInterval (index 0 = before first event -> N = 0)
if (nrow(n_tbl) == 0) {
  df$N <- 0
} else {
  # first Fertiliser == 1 half-hour per event date, or midnight if no dummy
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
  
  idx <- findInterval(df$timestamp, n_starts)
  df$N <- if_else(idx == 0L, 0, as.numeric(n_vals[pmax(idx, 1L)]))
}

# save updated RDS in place
saveRDS(df, in_path)
