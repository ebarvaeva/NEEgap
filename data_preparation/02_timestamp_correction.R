# =============================================================================
# 02_timestamp_correction.R — Timestamp Standardisation and Full-Year Grid
# =============================================================================
#
# PURPOSE
#   Takes the QC-passed merged CSV and imposes a complete, gapless half-hourly
#   timestamp grid: YYYY-01-01 00:00 → YYYY-12-31 23:30 (UTC, 30 min).
#   Half-hours absent from the QC file become NA rows.  Duplicate timestamps
#   are removed, keeping the first occurrence.
#
# WHY THIS STEP IS NECESSARY
#   Downstream scripts require a regular, gapless time axis.  Without it,
#   row-position indexing is unreliable and rolling-window computations
#   (rain_rolling_24, Phytomass Index) produce incorrect values near gaps.
#
# INPUTS
#   data/data_qc/{SITE}_{YEAR}/merged_qc.csv
#     Must contain TIMESTAMP_START (FluxNet integer YYYYMMDDHHMM)
#     or a 'timestamp' column parseable as POSIXct.
#
# OUTPUTS
#   data/data_qc/{SITE}_{YEAR}/merged_qc_fullgrid.csv
#     Complete 30-min grid; rows with no measurement have NA for all data columns.
#
# SECTIONS
#   User settings
#   1  Load QC file      -- parse timestamp; handle FluxNet and standard formats
#   2  Deduplicate       -- remove duplicate timestamps (keep first occurrence)
#   3  Full-year grid    -- build reference sequence of 17 520 half-hours
#   4  Left-join         -- merge QC data onto grid; gaps become NA rows
#   5  Save output
#
# NOTE ON TIMESTAMP CONVENTION
#   TIMESTAMP_START marks the START of the 30-min averaging period.
#   Row 1 (00:00) covers 00:00-00:30; row 17520 (23:30) covers 23:30-00:00.
#   All downstream scripts use this convention.
#
# =============================================================================

library(here)
library(dplyr)
library(readr)
library(lubridate)
library(tidyr)

# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------
# When sourced by 01_run_data_preparation.R, SITE and YEAR are injected
# automatically.  The defaults below allow standalone use during development.
if (!exists("SITE")) SITE <- "JC1"   # change for standalone use
if (!exists("YEAR")) YEAR <- 2023    # change for standalone use
if (!exists("TZ"))   TZ   <- "UTC"
TZ   <- "UTC"   # all timestamps stored and processed in UTC throughout

IN_DIR  <- here::here("data", "data_qc", paste0("results_qc_", SITE, "_", YEAR))
OUT_DIR <- IN_DIR   # write into the same folder

# -----------------------------------------------------------------------------
# SECTION 1 — Load QC file and parse timestamp
# -----------------------------------------------------------------------------
# Two timestamp formats are supported:
#   TIMESTAMP_START -- FluxNet integer format YYYYMMDDHHMM, converted to POSIXct UTC.
#   timestamp       -- POSIXct string written by Script 01 in Europe/Dublin time;
#                      converted to UTC for alignment with the reference grid.
# Rows outside the target year are discarded after parsing (can arise from
# EddyPro output files that bleed across the year boundary).
in_path <- file.path(IN_DIR, "merged_qc.csv")
if (!file.exists(in_path)) stop("QC file not found: ", in_path)

message("Loading ", in_path, " ...")
df_qc <- read_csv(in_path, show_col_types = FALSE)

# Parse timestamp ─────────────────────────────────────────────────────────────
# FluxNet TIMESTAMP_START: integer format YYYYMMDDHHMM
if ("TIMESTAMP_START" %in% names(df_qc) && !"timestamp" %in% names(df_qc)) {
  df_qc <- df_qc %>%
    mutate(timestamp = as.POSIXct(
      as.character(TIMESTAMP_START), format = "%Y%m%d%H%M", tz = TZ
    ))
} else if ("timestamp" %in% names(df_qc)) {
  df_qc <- df_qc %>%
    mutate(
      # Script 01 creates timestamps in Europe/Dublin; read them as such
      # then convert to UTC so the full-year UTC grid join aligns correctly
      timestamp = with_tz(
        as.POSIXct(timestamp, tz = "Europe/Dublin"),
        "UTC"
      )
    )
} else {
  stop("No usable timestamp column found (need TIMESTAMP_START or timestamp).")
}

# Filter to the target year only (removes any bleed from previous/next year)
df_qc <- df_qc %>%
  filter(year(timestamp) == YEAR) %>%
  arrange(timestamp)

message("  Rows in QC file for ", YEAR, ": ", nrow(df_qc))

# -----------------------------------------------------------------------------
# SECTION 2 — Remove duplicate timestamps
# -----------------------------------------------------------------------------
# Duplicates arise from overlapping EddyPro output files after restarts.
# The first occurrence is retained on the assumption it represents the
# original, uninterrupted measurement run.
n_before <- nrow(df_qc)
df_qc <- df_qc %>%
  distinct(timestamp, .keep_all = TRUE)
n_removed <- n_before - nrow(df_qc)
if (n_removed > 0) {
  message("  Removed ", n_removed, " duplicate timestamp(s) (kept first occurrence).")
}

# -----------------------------------------------------------------------------
# SECTION 3 — Build full-year 30-min reference grid
# -----------------------------------------------------------------------------
# The reference grid is a regular sequence of 17 520 half-hours (17 568 in
# a leap year) from 00:00 on 1 January to 23:30 on 31 December.
# A stopifnot() assertion verifies the count before the join.
start_ts <- as.POSIXct(paste0(YEAR, "-01-01 00:00:00"), tz = TZ)
end_ts   <- as.POSIXct(paste0(YEAR, "-12-31 23:30:00"), tz = TZ)
full_grid <- tibble(
  timestamp = seq(from = start_ts, to = end_ts, by = "30 min")
)

expected_rows <- as.integer(difftime(end_ts, start_ts, units = "mins") / 30) + 1L
stopifnot(nrow(full_grid) == expected_rows)
message("  Expected half-hours in ", YEAR, ": ", nrow(full_grid),
        " (", expected_rows, ")")

# -----------------------------------------------------------------------------
# SECTION 4 — Left-join QC data onto the full-year grid
# -----------------------------------------------------------------------------
# A left join with the reference grid as the left table inserts NA rows for
# every half-hour absent from the QC file.  These NA rows represent genuine
# measurement gaps; the gap-filling models will predict values for them later.
# The stopifnot() at the end guards against accidental row duplication.
df_full <- full_grid %>%
  left_join(df_qc, by = "timestamp") %>%
  arrange(timestamp)

n_missing <- sum(is.na(df_full$co2_flux_base_filters_6_70_grid) &
                   !df_full$timestamp %in% df_qc$timestamp)
message("  Half-hours added from grid (instrument gaps): ",
        nrow(df_full) - nrow(df_qc))
message("  Total rows in output: ", nrow(df_full))

# Quick sanity check: no duplicate timestamps in output
stopifnot(nrow(df_full) == nrow(distinct(df_full, timestamp)))

# -----------------------------------------------------------------------------
# SECTION 5 — Save output
# -----------------------------------------------------------------------------
# Written to the same folder as the input QC file so that Script 03 can
# locate it via the standard path convention.
out_path <- file.path(OUT_DIR, "merged_qc_fullgrid.csv")
write_csv(df_full, out_path)
message("Saved: ", out_path)
message("\nScript 02 complete.")
