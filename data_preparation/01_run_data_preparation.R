# =============================================================================
# 01_run_data_preparation.R — Data Preparation Master Script
# =============================================================================
# PURPOSE
#   Orchestrates the complete data preparation pipeline for both sites (JC1
#   and JC2) and both years (2023, 2024), producing the final analysis-ready
#   datasets data/data_prepared/JC1.rds and data/data_prepared/JC2.rds.
#
# PIPELINE OVERVIEW
#   For each site × year combination (JC1: 2023, 2024; JC2: 2023, 2024):
#
#   02_timestamp_correction.R
#       Impose complete 30-min calendar-year grid; add NA rows for gaps.
#       Output → results/qc/{SITE}_{YEAR}/merged_qc_fullgrid.csv
#
#   03_predictors_meteo.R
#       PPFD (biomet + fallbacks + RF gap-fill), Rg, Temp, RH, rain, VPD,
#       temporal encodings, season dummies, night flag, rain_rolling_24.
#       Output → results/prepared/{SITE}_{YEAR}_predictors.rds
#
#   04_management_days_since.R
#       Grazing_days_since, Fertiliser_days_since.
#       Updates → results/prepared/{SITE}_{YEAR}_predictors.rds
#
#   05_nitrogen.R
#       N (kg N ha⁻¹) step function from fertiliser event schedule.
#       Updates → results/prepared/{SITE}_{YEAR}_predictors.rds
#
#   06_grass_canopy.R
#       grass_height, grass_biomass from sparse field measurements.
#       Updates → results/prepared/{SITE}_{YEAR}_predictors.rds
#
#   07_bind_years.R
#       Concatenate 2023 + 2024, validate, save final site RDS.
#       Output → data/data_prepared/{SITE}.rds
#
# PREREQUISITE
#   Script 01 (quality_control_nee/) must be run manually for each site × year
#   before this pipeline.  It produces data/data_qc/results_qc_{SITE}_{YEAR}/
#   merged_qc.csv, which is the input to Script 02.
#   Any site × year combination for which this file is absent is automatically
#   skipped with an informative message.
#
# HOW TO RUN
#   Open this script in RStudio with the project root as the working directory.
#   Source the file:  source("data_preparation/01_run_data_preparation.R")
#
# REQUIREMENTS
#   library(here)  — all sub-scripts use here::here() for path resolution.
#   All raw data files must be in the locations described in DATA LAYOUT.
#
# DATA LAYOUT (files you supply; not included in the repository)
# =============================================================================
# data/
#   raw_data/
#     {SITE}_{YEAR}_eddypro.csv   — EddyPro full-output CSV
#     {SITE}_{YEAR}_biomet.csv    — Biomet half-hourly CSV
#     {SITE}_{YEAR}_fluxnet.csv   — FluxNet-format CSV
#     {SITE}_{YEAR}_meta.csv      — Site metadata / config
#
#   site_polygon/
#     {SITE}_polygon.csv          — Site boundary coordinates
#                                   Columns: lon, lat (closed ring)
#
#   meteireann_data/
#     met_hourly.csv              — Hourly weather from Met Éireann
#                                   Columns: date (DD/MM/YYYY HH:MM), temp, rhum, rain
#     solar_hourly.xlsx           — Hourly global radiation
#                                   Columns: date, glorad (W m⁻²)
#
#   alt_ppfd/
#     alt_ppfd_{SITE}_{YEAR}.xlsx — Alternative PPFD from nearby stations
#                                   Sheet 1: Oakpark  (date, time, PPFD_1_1_1)
#                                   Sheet 2: {SITE}   (date, time, PPFD_1_1_1)
#                                   Sheet 3: CD_grass (date, time, PPFD_1_1_1)
#
#   management_data/
#     grazing_events_{SITE}.csv   — Grazing event windows
#                                   Columns: date, start_time, end_time, event_id
#     fertiliser_events_{SITE}.csv — Fertiliser event windows
#                                   Columns: date, start_time, end_time, event_id
#     last_events_before_year_{SITE}.csv
#                                   Columns: year_for, event_type, last_end
#                                   (one row per event_type per year to process)
#     nitrogen_amounts_{SITE}.csv — N kg ha⁻¹ per fertiliser event
#                                   Columns: fert_date (DD/MM/YYYY),
#                                            N_kg_ha, is_slurry (TRUE/FALSE)
#
#   canopy_data/
#     canopy_measurements_{SITE}_{YEAR}.csv
#                                   Columns: date, subplot_id,
#                                            height_cm, drymass_kg_ha
#
#   data_prepared/               — created automatically by Script 07
# =============================================================================

library(here)
source("data_preparation/generate_management_csvs.R")
# =============================================================================
# CONFIGURATION
# =============================================================================
# Edit these settings for a different site set, year range, or project layout.

SITES <- c("JC1", "JC2")   # sites to process; add more as needed
YEARS <- c(2023, 2024)      # years to process per site

# Set to FALSE for any site without management data (grazing, fertilisation).
# When FALSE, Scripts 04, 05, and 06 are skipped for that site and the final
# dataset will not contain management predictor columns.
HAS_MANAGEMENT <- list(JC1 = TRUE, JC2 = TRUE)

# Script paths
SCRIPT_DIR <- here::here("data_preparation")

scripts <- list(
  `02` = file.path(SCRIPT_DIR, "02_timestamp_correction.R"),
  `03` = file.path(SCRIPT_DIR, "03_predictors_meteo.R"),
  `04` = file.path(SCRIPT_DIR, "04_management_days_since.R"),
  `05` = file.path(SCRIPT_DIR, "05_nitrogen.R"),
  `06` = file.path(SCRIPT_DIR, "06_grass_canopy.R"),
  `07` = file.path(SCRIPT_DIR, "07_bind_years.R")
)

# Verify all scripts exist
missing_scripts <- names(scripts)[!vapply(scripts, file.exists, logical(1))]
if (length(missing_scripts))
  
  stop("Scripts not found: ", paste(missing_scripts, collapse=", "))

# =============================================================================
# HELPER: run_script()
# =============================================================================
# Sources a sub-script in an isolated environment, injecting SITE, YEAR,
# TWIN_SITE, and HAS_MANAGEMENT as variables.  Executing in a new environment
# prevents sub-script variables from accumulating in the global workspace and
# allows multiple site × year combinations to be processed in a single session
# without variable conflicts.
# twin_site is the other site's identifier, passed to Script 03 as the PPFD
# fallback source.

run_script <- function(script_path, site, year, twin_site = NULL, has_mgmt = TRUE) {
  env <- new.env(parent = globalenv())
  env$SITE      <- site
  env$YEAR      <- year
  if (!is.null(twin_site)) env$TWIN_SITE <- twin_site
  env$HAS_MANAGEMENT <- has_mgmt
  message("\n", strrep("─", 60))
  message("Running: ", basename(script_path))
  message("  SITE = ", site, "  |  YEAR = ", year)
  message(strrep("─", 60))
  source(script_path, local = env)
  invisible(NULL)
}

# =============================================================================
# MAIN PIPELINE
# =============================================================================
# The outer loop processes one site at a time.  The inner loop processes each
# year for that site in sequence (Scripts 02–06), then Script 07 binds the
# per-year outputs into the final site-level dataset.
#
# Error handling: if any script fails for a given year, the remaining scripts
# for that year are skipped (tryCatch), but processing continues with the next
# year.  Script 07 is called even if some years failed, and it will bind only
# the years whose RDS files are present.

for (site in SITES) {
  twin <- setdiff(SITES, site)   # the other site; used as PPFD fallback in Script 03
  mgmt <- isTRUE(HAS_MANAGEMENT[[site]])

  for (yr in YEARS) {
    # Check that the QC output exists before running Scripts 02–06.
    # If merged_qc.csv is absent (QC not yet run for this year), the entire
    # year is skipped.  Run quality_control_nee/ manually first.
    qc_input <- here::here("data", "data_qc",
                            paste0("results_qc_", site, "_", yr), "merged_qc.csv")
    if (!file.exists(qc_input)) {
      message("\n  QC output not found for ", site, " ", yr, " — skipping.")
      message("  Expected: ", qc_input)
      next
    }

    tryCatch({
      run_script(scripts[["02"]], site = site, year = yr)
      run_script(scripts[["03"]], site = site, year = yr, twin_site = twin)

      if (mgmt) {
        run_script(scripts[["04"]], site = site, year = yr)
        run_script(scripts[["05"]], site = site, year = yr)
        run_script(scripts[["06"]], site = site, year = yr)
      } else {
        message("Skipping management scripts (HAS_MANAGEMENT = FALSE for ", site, ")")
      }
    }, error = function(e) {
      message("\n  ERROR in ", site, " ", yr, ": ", conditionMessage(e))
      message("  Skipping remaining scripts for this year.")
    })
  }

  # --- Bind 2023 + 2024 for this site ---
  run_script(scripts[["07"]], site = site, year = NA, has_mgmt = mgmt)
}

# =============================================================================
# SUMMARY
# =============================================================================
# Prints the row count, column count, and number of finite NEE observations
# for each output dataset.  A missing output file indicates that Script 07
# failed or was not reached for that site — check the error messages above.

message("\n", strrep("=", 60))
message("DATA PREPARATION COMPLETE")
message(strrep("=", 60))

for (site in SITES) {
  out <- here::here("data", "data_prepared", paste0(site, ".rds"))
  if (file.exists(out)) {
    df_check <- readRDS(out)
    message(sprintf("  %-5s: %6d rows  ×  %d cols  |  finite NEE: %d",
                    site, nrow(df_check), ncol(df_check),
                    sum(is.finite(df_check$NEE_orig))))
    rm(df_check)
  } else {
    message("  ", site, ": OUTPUT FILE MISSING — check errors above.")
  }
}
