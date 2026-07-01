# 02_timestamp_correction.R — build a complete 30-min grid for one site x year
#
# Loads the QC output for SITE x YEAR, parses timestamps, drops duplicates, and
# left-joins onto a full calendar-year half-hourly grid so instrument gaps
# become explicit NA rows. JC1 2020 is stored as merged_qc.rds; every other
# site-year is read from merged_qc.csv. Output is always merged_qc_fullgrid.csv.
#
# Input  : data/data_qc/results_qc_{SITE}_{YEAR}/merged_qc.{csv,rds}
# Output : data/data_qc/results_qc_{SITE}_{YEAR}/merged_qc_fullgrid.csv

library(here)
library(dplyr)
library(readr)
library(lubridate)
library(tidyr)

# settings (SITE, YEAR injected by script 01)
if (!exists("SITE")) SITE <- "JC1"
if (!exists("YEAR")) YEAR <- 2023
TZ   <- "UTC"

IN_DIR  <- here::here("data", "data_qc", paste0("results_qc_", SITE, "_", YEAR))
OUT_DIR <- IN_DIR

# load QC file: JC1 2020 is an .rds, every other site-year is a .csv
is_2020_jc1 <- (SITE == "JC1" && YEAR == 2020)
if (is_2020_jc1) {
  in_path <- file.path(IN_DIR, "merged_qc.rds")
  if (!file.exists(in_path)) stop("QC RDS not found: ", in_path)
  df_qc <- readRDS(in_path)
} else {
  in_path <- file.path(IN_DIR, "merged_qc.csv")
  if (!file.exists(in_path)) stop("QC CSV not found: ", in_path)
  df_qc <- read_csv(in_path, show_col_types = FALSE)
}

# build a POSIXct timestamp from whichever column is available
if ("TIMESTAMP_START" %in% names(df_qc) && !"timestamp" %in% names(df_qc)) {
  df_qc <- df_qc %>%
    mutate(timestamp = as.POSIXct(
      as.character(TIMESTAMP_START), format = "%Y%m%d%H%M", tz = TZ
    ))
} else if ("timestamp" %in% names(df_qc)) {
  df_qc <- df_qc %>%
    mutate(timestamp = with_tz(as.POSIXct(timestamp, tz = "UTC"), "UTC"))
} else {
  stop("No usable timestamp column found (need TIMESTAMP_START or timestamp).")
}

# keep only rows in the target year, ordered by time
df_qc <- df_qc %>%
  filter(year(timestamp) == YEAR) %>%
  arrange(timestamp)

# drop duplicate timestamps, keeping the first occurrence
df_qc <- df_qc %>% distinct(timestamp, .keep_all = TRUE)

# full-year 30-min grid from Jan 1 00:00 to Dec 31 23:30
start_ts <- as.POSIXct(paste0(YEAR, "-01-01 00:00:00"), tz = TZ)
end_ts   <- as.POSIXct(paste0(YEAR, "-12-31 23:30:00"), tz = TZ)
full_grid <- tibble(timestamp = seq(from = start_ts, to = end_ts, by = "30 min"))

expected_rows <- as.integer(difftime(end_ts, start_ts, units = "mins") / 30) + 1L
stopifnot(nrow(full_grid) == expected_rows)

# left-join QC data onto the grid so instrument gaps become NA rows
df_full <- full_grid %>%
  left_join(df_qc, by = "timestamp") %>%
  arrange(timestamp)

stopifnot(nrow(df_full) == nrow(distinct(df_full, timestamp)))

# save the full-grid CSV (downstream scripts always read this)
out_path <- file.path(OUT_DIR, "merged_qc_fullgrid.csv")
write_csv(df_full, out_path)
