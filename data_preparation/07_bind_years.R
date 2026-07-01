# 07_bind_years.R — bind years and save the final per-site RDS
#
# Loads each available year's predictor RDS for SITE, keeps common columns,
# binds them into one time-ordered data frame, validates it (no duplicate
# timestamps, expected full-year row count, zero NA in required predictors),
# and saves the final {SITE}.rds. JC1 binds 2020 + 2023 + 2024; JC2 binds
# 2023 + 2024.
#
# Input  : data/data_prepared/{SITE}_{YEAR}_predictors.rds  (one per available year)
# Output : data/data_prepared/{SITE}.rds

library(here)
library(dplyr)
library(lubridate)

# settings (SITE injected by script 01; YEAR is NA for this script)
if (!exists("SITE")) SITE <- "JC1"

# site-specific year list — JC1 also has 2020
YEARS <- if (SITE == "JC1") c(2020, 2023, 2024) else c(2023, 2024)

REQUIRED_PREDICTORS <- c(
  "PPFD", "Rg", "VPD", "RH", "Temp", "rain", "rain_rolling_24", "night",
  "hour_sin", "hour_cos", "doy_sin", "doy_cos", "month_sin", "month_cos",
  "Winter", "Spring", "Summer", "Autumn"
)

MANAGEMENT_PREDICTORS <- c(
  "Grazing_days_since", "Fertiliser_days_since", "N",
  "grass_height", "grass_biomass"
)

if (!exists("HAS_MANAGEMENT")) HAS_MANAGEMENT <- TRUE

OUT_DIR <- here::here("data/data_prepared")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# load each available year (missing years are skipped)
dfs <- lapply(YEARS, function(yr) {
  p <- here::here("data", "data_prepared",
                  paste0(SITE, "_", yr, "_predictors.rds"))
  if (!file.exists(p)) return(NULL)
  df_yr <- readRDS(p) %>%
    mutate(timestamp = as.POSIXct(timestamp, tz = "UTC")) %>%
    arrange(timestamp) %>%
    distinct(timestamp, .keep_all = TRUE)
  df_yr
})

dfs <- Filter(Negate(is.null), dfs)
if (length(dfs) == 0)
  stop("No predictor RDS files found for ", SITE,
       ". Run QC + pipeline for at least one year.")

# keep only columns common to all years, then bind
common_cols <- Reduce(intersect, lapply(dfs, names))
dfs <- lapply(dfs, function(d) d[, common_cols, drop = FALSE])

df <- dplyr::bind_rows(dfs) %>%
  arrange(timestamp) %>%
  distinct(timestamp, .keep_all = TRUE)

# no duplicate timestamps after bind
n_dup <- sum(duplicated(df$timestamp))
if (n_dup > 0)
  stop(n_dup, " duplicate timestamps remain after bind.")

# row count should match a full 30-min grid summed over the available years
available_years <- unique(year(df$timestamp))
expected <- sum(vapply(available_years, function(yr) {
  start_ts <- as.POSIXct(paste0(yr, "-01-01 00:00:00"), tz = "UTC")
  end_ts   <- as.POSIXct(paste0(yr, "-12-31 23:30:00"), tz = "UTC")
  as.integer(difftime(end_ts, start_ts, units = "mins") / 30) + 1L
}, integer(1)))

if (nrow(df) != expected)
  warning("Row count differs from expected full-year grid by ",
          nrow(df) - expected, " rows. Check for gaps or extra rows.")

# required columns must all be present
required_cols <- c("timestamp", "NEE_orig", REQUIRED_PREDICTORS)
if (HAS_MANAGEMENT) required_cols <- c(required_cols, MANAGEMENT_PREDICTORS)

missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols))
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

# required predictors (excluding optional canopy columns) must contain no NA
OPTIONAL_COLS <- c("grass_height", "grass_biomass")
check_cols    <- setdiff(required_cols, c("timestamp", "NEE_orig", OPTIONAL_COLS))
na_tbl <- tibble::tibble(
  column = check_cols,
  n_na   = vapply(check_cols, function(cc) sum(is.na(df[[cc]])), integer(1))
) %>% dplyr::arrange(dplyr::desc(n_na))

any_na <- any(na_tbl$n_na > 0)
if (any_na)
  stop("Some predictors still contain NA values. ",
       "Re-run the relevant preparation scripts or investigate missing data.")

# select final columns and save
df_final <- df %>%
  select(all_of(required_cols)) %>%
  arrange(timestamp)

out_path <- file.path(OUT_DIR, paste0(SITE, ".rds"))
saveRDS(df_final, out_path)
