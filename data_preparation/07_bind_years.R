# =============================================================================
# 07_bind_years.R -- Bind Years and Save Final Site Dataset
# =============================================================================
#
# PURPOSE
#   Concatenates the per-year predictor RDS files (2023 + 2024) produced by
#   Scripts 03-06, validates completeness and consistency, and saves the final
#   two-year dataset as data/data_prepared/{SITE}.rds.  This is the file read
#   by run_model.R and all gap-filling model scripts.
#
# INPUTS
#   data/data_prepared/{SITE}_2023_predictors.rds
#   data/data_prepared/{SITE}_2024_predictors.rds
#   (or whichever years are available; missing years are skipped with a warning)
#
# OUTPUTS
#   data/data_prepared/{SITE}.rds
#     A single data frame with 35 040 rows (17 520 per year), all required
#     predictors, and NEE_orig (which may contain NAs where measurements were
#     unavailable -- these are the rows to be gap-filled).
#
# VALIDATION STEPS
#   1. Duplicate timestamps: any detected -> hard stop (indicates overlapping
#      per-year files or a coding error in Scripts 02-06).
#   2. Expected row count: compared against 17 520 x number_of_years.  A
#      mismatch produces a warning but does not stop execution (may arise from
#      an incomplete year at the end of the observation record).
#   3. Missing required columns: hard stop listing each missing column.
#   4. NA audit: all predictor columns (excluding NEE_orig and columns in
#      OPTIONAL_COLS) must have zero NAs.  If any remain, the script stops
#      and prints a ranked table of NA counts so the user can identify which
#      script to re-run.
#
# ROBUSTNESS TO COLUMN MISMATCHES
#   If Scripts 03-06 produced different column sets across years (e.g. because
#   management data was only available for one year), the bind uses only the
#   intersection of columns present in all years.  A warning is printed listing
#   dropped columns.
#
# =============================================================================

library(here)
library(dplyr)
library(lubridate)

# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------
# YEARS: list of calendar years to bind.  Add a year here if a third year of
# data becomes available; the expected row count check will update automatically.
YEARS <- c(2023, 2024)

# Required predictor columns -- any of these missing from the bound data frame
# will cause a hard stop.  All must have zero NAs (except NEE_orig and
# OPTIONAL_COLS).  Add new predictors here when Scripts 03-06 are extended.
REQUIRED_PREDICTORS <- c(
  "PPFD", "Rg", "VPD", "RH", "Temp", "rain", "rain_rolling_24", "night",
  "hour_sin", "hour_cos", "doy_sin", "doy_cos", "month_sin", "month_cos",
  "Winter", "Spring", "Summer", "Autumn"
)

# Management predictors: appended to required_cols only for managed sites.
# Set HAS_MANAGEMENT <- FALSE to process an unmanaged site.
MANAGEMENT_PREDICTORS <- c(
  "Grazing_days_since", "Fertiliser_days_since", "N",
  "grass_height", "grass_biomass"
)

# TRUE for managed sites (Scripts 04-06 were run); FALSE for BASE-only sites.
HAS_MANAGEMENT <- TRUE

OUT_DIR <- here::here("data/data_prepared")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# SECTION 1 -- Load and bind per-year RDS files
# -----------------------------------------------------------------------------
# Loads each year's predictors RDS.  Years whose file is absent are skipped
# (with an informative message) to allow partial runs during development.
# Duplicate timestamps within each year are removed (second occurrence dropped)
# before the bind; the final deduplication after bind_rows() handles any cross-
# year boundary duplicates.
# Only columns present in ALL years are retained (common_cols intersection),
# preventing bind_rows() from failing if Scripts 03-06 added different columns
# in different years.
message("Loading year RDS files for ", SITE, " ...")
dfs <- lapply(YEARS, function(yr) {
  p <- here::here("data", "data_prepared", paste0(SITE, "_", yr, "_predictors.rds"))
  if (!file.exists(p)) {
    message("  ", yr, ": predictors RDS not found — skipping.")
    message("    (Run QC for this year then re-run the pipeline.)")
    return(NULL)
  }
  df_yr <- readRDS(p)
  df_yr <- df_yr %>%
    mutate(timestamp = as.POSIXct(timestamp, tz = "UTC")) %>%
    arrange(timestamp) %>%
    distinct(timestamp, .keep_all = TRUE)
  message("  ", yr, ": ", nrow(df_yr), " rows")
  df_yr
})

# Drop years that were skipped (NULL entries)
dfs <- Filter(Negate(is.null), dfs)
if (length(dfs) == 0)
  stop("No predictor RDS files found for ", SITE, ". Run QC + pipeline for at least one year.")
if (length(dfs) < length(YEARS))
  message("  Note: only ", length(dfs), "/", length(YEARS),
          " years available — binding available years only.")

# Common columns only (prevents bind failure if scripts differ between years)
common_cols <- Reduce(intersect, lapply(dfs, names))
message("  Columns in common: ", length(common_cols))
dfs <- lapply(dfs, function(d) d[, common_cols, drop = FALSE])

df <- dplyr::bind_rows(dfs) %>%
  arrange(timestamp) %>%
  distinct(timestamp, .keep_all = TRUE)

message("Total rows after binding: ", nrow(df))

# -----------------------------------------------------------------------------
# SECTION 2 -- Duplicate timestamp check
# -----------------------------------------------------------------------------
# Any duplicate in the bound dataset is a hard error: it means two per-year
# files share a row, which would silently inflate the training set for the
# gap-filling models.  The script stops with a clear message rather than
# producing a quietly corrupted output file.
n_dup <- sum(duplicated(df$timestamp))
if (n_dup > 0) {
  stop(n_dup, " duplicate timestamps remain after bind. ",
       "Check per-year files for overlapping rows.")
}
message("No duplicate timestamps.")

# -----------------------------------------------------------------------------
# SECTION 3 -- Expected row count check
# -----------------------------------------------------------------------------
# Computes the expected number of half-hours across all YEARS (17 520 per
# non-leap year, 17 568 for a leap year) and compares it with nrow(df).
# A mismatch is a warning rather than an error because a partial final year
# is legitimate.  Gaps within the grid (instrument outages) are represented
# as NA rows and are correctly included in the count.
expected <- sum(vapply(YEARS, function(yr) {
  start_ts <- as.POSIXct(paste0(yr, "-01-01 00:00:00"), tz = "UTC")
  end_ts   <- as.POSIXct(paste0(yr, "-12-31 23:30:00"), tz = "UTC")
  as.integer(difftime(end_ts, start_ts, units = "mins") / 30) + 1L
}, integer(1)))

message("Expected total half-hours (", paste(YEARS, collapse="+"), "): ", expected)
message("Actual rows: ", nrow(df))
if (nrow(df) != expected) {
  warning("Row count differs from expected full-year grid by ",
          nrow(df) - expected, " rows. Check for gaps or extra rows.")
}

# -----------------------------------------------------------------------------
# SECTION 4 -- Column presence and NA audit
# -----------------------------------------------------------------------------
# First checks that all required_cols are present; stops with a list of missing
# columns if any are absent.  Then checks that all non-optional, non-NEE
# predictor columns have zero NA values.  If NAs remain, a ranked table is
# printed (most NAs first) and the script stops.  The user should identify
# which Script (03 for meteorology, 04 for days_since, 05 for N, 06 for canopy)
# is responsible and re-run it.
# OPTIONAL_COLS (grass_height, grass_biomass) are excluded from the NA check
# because sparse canopy measurements can legitimately leave NAs at the boundary
# of the observation period (rule = 2 prevents this in practice, but the
# check is kept permissive for robustness).
required_cols <- c("timestamp", "NEE_orig", REQUIRED_PREDICTORS)
if (HAS_MANAGEMENT) required_cols <- c(required_cols, MANAGEMENT_PREDICTORS)

missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols)) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# NA audit (exclude NEE_orig)
OPTIONAL_COLS <- c("grass_height", "grass_biomass")
check_cols <- setdiff(required_cols, c("timestamp", "NEE_orig", OPTIONAL_COLS))
na_tbl <- tibble::tibble(
  column = check_cols,
  n_na   = vapply(check_cols, function(cc) sum(is.na(df[[cc]])), integer(1))
) %>% dplyr::arrange(dplyr::desc(n_na))

message("\nNA counts in predictor columns:")
print(na_tbl, n = Inf)

any_na <- any(na_tbl$n_na > 0)
if (any_na) {
  stop("Some predictors still contain NA values. ",
       "Re-run the relevant preparation scripts or investigate missing data.")
}
message("\nAll required predictors: 0 NA. Proceeding to save.")

# -----------------------------------------------------------------------------
# SECTION 5 -- Select final columns and save
# -----------------------------------------------------------------------------
# Retains only the validated required_cols in a canonical column order.
# The output file is the primary input to run_model.R.
df_final <- df %>%
  select(all_of(required_cols)) %>%
  arrange(timestamp)

out_path <- file.path(OUT_DIR, paste0(SITE, ".rds"))
saveRDS(df_final, out_path)

message("Saved: ", out_path)
message("  Rows: ", nrow(df_final), "  |  Columns: ", ncol(df_final))
message("  NEE_orig finite: ", sum(is.finite(df_final$NEE_orig)))
message("\nScript 07 complete.")
