# =============================================================================
# Regrowth_period.R — Regrowth Period and Animal/Nitrogen Dummies (SKELETON)
# =============================================================================
#
# PURPOSE
#   Uses the management event tables from ManagementEvents.R to add contextual
#   columns to the gap-filling data frame:
#
#   grazing_group_{i}   — 1 while the i-th grazing bout is in progress
#   growth_{i}          — 1 during the regrowth interval after grazing bout i
#   pre_grazing         — 1 before the first grazing event of the record
#   post_grazing        — 1 after the last grazing event of the record
#   cows                — number of animals present (0 outside grazing windows)
#   N                   — kg N ha⁻¹ applied at the most recent fertiliser event
#                         (step function, 0 before first application)
#
# HOW TO FILL IN YOUR DATA
#   1. Fill in cow_windows with the start date, end date, and animal count for
#      each grazing period.  One row per continuous occupancy window.
#   2. Fill in N_windows with the date and N amount for each mineral fertiliser
#      application.  Slurry events should be excluded (see comment in table).
#   3. Update the three "CHANGE SITE" lines in Section 2 to point to the correct
#      site identifier, RDS file, and event tables from ManagementEvents.R.
#
# PLACEHOLDER CONVENTION
#   "YYYY-MM-DD"  — replace with an actual date string
#   0L            — replace with the actual animal count (integer)
#   "DD/MM/YYYY"  — replace with a date in day/month/year format
#   0.0           — replace with the actual N amount (kg N ha⁻¹)
#
# SOURCED AFTER
#   ManagementEvents.R  (must be sourced first; provides the *_events_* tables)
#
# =============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(here)

# ManagementEvents.R must already be sourced (or sourced here):
source(here::here("data_preparation", "ManagementEvents.R"))


# =============================================================================
# SECTION 1 — Site-level data tables  (FILL IN YOUR DATA HERE)
# =============================================================================

# ---------------------------------------------------------------------------
# cow_windows
# ---------------------------------------------------------------------------
# One row per continuous animal-occupancy window.  Windows should align with
# the grazing events in ManagementEvents.R — use the same start and end dates.
#
# Columns:
#   start_date  — "YYYY-MM-DD"  first day animals were present
#   end_date    — "YYYY-MM-DD"  last day animals were present (inclusive)
#   cow_number  — integer       number of animals in this window

cow_windows_site1 <- tibble::tribble(
  ~start_date,    ~end_date,      ~cow_number,
  # "YYYY-MM-DD",  "YYYY-MM-DD",  0L,   # <-- replace with your first grazing window
  # "YYYY-MM-DD",  "YYYY-MM-DD",  0L,   # <-- add one row per grazing bout
  "YYYY-MM-DD",  "YYYY-MM-DD",   0L    # placeholder — remove and replace
)

cow_windows_site2 <- tibble::tribble(
  ~start_date,    ~end_date,      ~cow_number,
  "YYYY-MM-DD",  "YYYY-MM-DD",   0L    # placeholder — remove and replace
)

# ---------------------------------------------------------------------------
# N_windows
# ---------------------------------------------------------------------------
# One row per mineral fertiliser application.
# IMPORTANT: slurry events must be excluded (uncertain N content).
# If a date appears in ManagementEvents.R as "Slurry", do NOT list it here.
#
# Columns:
#   fert_date  — "DD/MM/YYYY"  date of application (day/month/year format)
#   N_kg_ha    — numeric       kg N ha⁻¹ applied

N_windows_site1 <- tibble::tribble(
  ~fert_date,     ~N_kg_ha,
  # "DD/MM/YYYY",  0.0,     # <-- replace with your first fertiliser event
  # "DD/MM/YYYY",  0.0,     # <-- add one row per mineral application
  # "DD/MM/YYYY"  # Slurry — excluded on purpose (do not include)
  "DD/MM/YYYY",   0.0      # placeholder — remove and replace
)

N_windows_site2 <- tibble::tribble(
  ~fert_date,     ~N_kg_ha,
  # "DD/MM/YYYY"  # Slurry — excluded on purpose
  "DD/MM/YYYY",   0.0      # placeholder — remove and replace
)


# =============================================================================
# SECTION 2 — Site selection  (CHANGE THESE THREE LINES PER SITE)
# =============================================================================
# Set site_id, the RDS file to augment, and which event and data tables to use.
# Run this script once per site (or wrap in a loop over SITES).

site_id     <- "SITE1"                 # <-- CHANGE SITE identifier
rds_name    <- "SITE1_2023_2024.rds"   # <-- CHANGE RDS file name

events_2023 <- site1_events_2023       # <-- CHANGE to matching events table from ManagementEvents.R
events_2024 <- site1_events_2024       # <-- CHANGE
cow_windows <- cow_windows_site1       # <-- CHANGE
N_windows   <- N_windows_site1         # <-- CHANGE


# =============================================================================
# SECTION 3 — Build grazing groups  (DO NOT EDIT)
# =============================================================================
# Groups consecutive grazing events that are separated by fewer than 10 days
# into a single "grazing bout".  Each bout gets a grazing_group index.
# The interval between bouts is the "growth" (regrowth) period.

grazing_events <- dplyr::bind_rows(events_2023, events_2024) %>%
  dplyr::filter(stringr::str_detect(event, "^Grazing")) %>%
  dplyr::mutate(
    date_dt  = as.POSIXct(date, tz = "UTC"),
    start_dt = date_dt +
      lubridate::hm(dplyr::if_else(start_time == "24:00", "00:00", start_time)) +
      lubridate::days(dplyr::if_else(start_time == "24:00", 1, 0)),
    end_dt = date_dt +
      lubridate::hm(dplyr::if_else(end_time == "24:00", "00:00", end_time)) +
      lubridate::days(dplyr::if_else(end_time == "24:00", 1, 0))
  ) %>%
  dplyr::arrange(start_dt) %>%
  dplyr::mutate(
    grazing_group = cumsum(
      dplyr::if_else(
        is.na(dplyr::lag(start_dt)) |
          as.numeric(difftime(start_dt, dplyr::lag(start_dt), units = "days")) > 10,
        1L, 0L
      )
    )
  )

grazing_groups <- grazing_events %>%
  dplyr::group_by(grazing_group) %>%
  dplyr::summarise(
    group_start    = min(start_dt),
    last_event_end = max(end_dt),
    .groups = "drop"
  ) %>%
  dplyr::arrange(group_start) %>%
  dplyr::mutate(group_end = dplyr::coalesce(dplyr::lead(group_start), last_event_end))


# =============================================================================
# SECTION 4 — Load data frame and add dummy columns  (DO NOT EDIT)
# =============================================================================
# Reads the prepared RDS for the selected site and adds all dummy columns.
# The RDS must already contain: timestamp, Grazing, Fertiliser (from scripts 03-04).

df <- readRDS(here::here("data", "data_prepared", site_id, rds_name))

tz_df <- attr(df$timestamp, "tzone")
if (is.null(tz_df) || tz_df == "") tz_df <- "UTC"

grazing_groups_tz <- grazing_groups %>%
  dplyr::mutate(
    grazing_group  = as.integer(grazing_group),
    group_start    = lubridate::with_tz(group_start,    tz_df),
    last_event_end = lubridate::with_tz(last_event_end, tz_df)
  ) %>%
  dplyr::arrange(group_start)

first_start <- grazing_groups_tz$group_start[1]
last_end    <- grazing_groups_tz$last_event_end[nrow(grazing_groups_tz)]

df <- df %>%
  dplyr::mutate(
    pre_grazing  = as.integer(timestamp < first_start),
    post_grazing = as.integer(timestamp >= last_end),
    .row_id      = dplyr::row_number()
  )

# --- grazing_group_{i}: 1 while bout i is in progress ----------------------
grazing_group_dummies <- df %>%
  dplyr::select(.row_id, timestamp) %>%
  tidyr::crossing(
    grazing_groups_tz %>% dplyr::select(grazing_group, group_start, last_event_end)
  ) %>%
  dplyr::mutate(
    val = as.integer(timestamp >= group_start & timestamp < last_event_end),
    var = paste0("grazing_group_", grazing_group)
  ) %>%
  dplyr::select(.row_id, var, val) %>%
  tidyr::pivot_wider(names_from = var, values_from = val, values_fill = 0)

# --- growth_{i}: 1 during regrowth interval after bout i ------------------
growth_def <- grazing_groups_tz %>%
  dplyr::transmute(
    grazing_group,
    growth_start = last_event_end,
    growth_end   = dplyr::lead(group_start)
  ) %>%
  dplyr::mutate(growth_end = dplyr::coalesce(growth_end, growth_start))

growth_dummies <- df %>%
  dplyr::select(.row_id, timestamp) %>%
  tidyr::crossing(growth_def) %>%
  dplyr::mutate(
    val = as.integer(timestamp >= growth_start & timestamp < growth_end),
    var = paste0("growth_", grazing_group)
  ) %>%
  dplyr::select(.row_id, var, val) %>%
  tidyr::pivot_wider(names_from = var, values_from = val, values_fill = 0)

df <- df %>%
  dplyr::left_join(grazing_group_dummies, by = ".row_id") %>%
  dplyr::left_join(growth_dummies,        by = ".row_id")

# --- cows: animal count during active grazing windows ---------------------
cow_windows2 <- cow_windows %>%
  dplyr::mutate(
    start_date = lubridate::ymd(start_date),
    end_date   = lubridate::ymd(end_date),
    start_dt   = as.POSIXct(start_date, tz = tz_df),
    end_dt     = as.POSIXct(end_date + lubridate::days(1), tz = tz_df)
  )

df <- df %>%
  dplyr::left_join(
    df %>%
      dplyr::select(.row_id, timestamp, Grazing) %>%
      tidyr::crossing(cow_windows2 %>% dplyr::select(start_dt, end_dt, cow_number)) %>%
      dplyr::mutate(hit = Grazing == 1 & timestamp >= start_dt & timestamp < end_dt) %>%
      dplyr::group_by(.row_id) %>%
      dplyr::summarise(
        cows = dplyr::if_else(any(hit), first(cow_number[hit]), 0L),
        .groups = "drop"
      ),
    by = ".row_id"
  )

# --- N: step function from mineral fertiliser applications ----------------
N_windows2 <- N_windows %>%
  dplyr::mutate(
    fert_date = lubridate::dmy(stringr::str_replace_all(fert_date, "\\.", "/"))
  )

# Pin the N step to the first observed half-hour where Fertiliser == 1 on
# each fertiliser date (fallback to midnight if the dummy is absent).
fert_starts <- df %>%
  dplyr::filter(Fertiliser == 1) %>%
  dplyr::mutate(fert_date = as.Date(timestamp)) %>%
  dplyr::group_by(fert_date) %>%
  dplyr::summarise(start_dt_obs = min(timestamp), .groups = "drop")

fert_schedule <- N_windows2 %>%
  dplyr::left_join(fert_starts, by = "fert_date") %>%
  dplyr::mutate(
    start_dt = dplyr::coalesce(start_dt_obs, as.POSIXct(fert_date, tz = tz_df))
  ) %>%
  dplyr::arrange(start_dt)

n_starts <- fert_schedule$start_dt
n_vals   <- fert_schedule$N_kg_ha

df <- df %>%
  dplyr::mutate(
    .n_idx  = findInterval(timestamp, n_starts),
    .n_idx2 = dplyr::if_else(.n_idx == 0L, 1L, .n_idx),
    N       = dplyr::if_else(.n_idx == 0L, 0, as.numeric(n_vals[.n_idx2]))
  ) %>%
  dplyr::select(-.n_idx, -.n_idx2)


# =============================================================================
# SECTION 5 — Save  (DO NOT EDIT)
# =============================================================================

out_path <- here::here("data", "data_prepared", site_id, rds_name)
saveRDS(df, out_path)
message("Saved: ", out_path)

# ====================== end Regrowth_period.R ================================
