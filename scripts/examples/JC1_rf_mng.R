# JC1_rf_mng.R — Run script for Random Forest NEE gap-filling
#
# Entry point for one site/management run. It loads the prepared dataset, builds
# the predictor list from FEATURE_SET, creates the output directory, and sources
# models/RF_CV.R to gap-fill NEE by cross-validation across four gap-size classes
# (VL = very large, L = large, M = medium, S = small).
#
# Required input columns (in DATA_FILE):
#   timestamp   — POSIXct / parseable half-hourly date-time
#   NEE_orig    — net ecosystem exchange (µmol m-2 s-1)
#   night       — binary flag (1 = night, 0 = day)
#   + any FEATURE_SET columns present (absent ones are silently dropped)
#   Grazing_days_since is optional; with it (and PI_ENABLED) the pre-computed
#   PI_{gap_label} columns from Script 09 are added as predictors.
#
# Outputs (in RESULTS_DIR):
#   df_cv_all_predictions.rds   — gap-filled NEE (4 prediction columns: VL/L/M/S)
#   variable_importance/
#     rf_vi_{SIZE}_{LABEL}.rds
#     rf_variable_importance_{SIZE}_all.{rds,csv}
#
# To adapt: change SITE_NAME / MANAGEMENT_CONDITION, point DATA_FILE at the new
# .rds, and edit FEATURE_SET as needed.


# Edit these settings for each site / run:
SITE_NAME            <- "JC1"        # site identifier ("JC2", "JC3", …)
MANAGEMENT_CONDITION <- "managed"    # "managed" or "unmanaged"
MODEL_CHOICE         <- "RF"         # names the model script (RF_CV.R) and results subfolder

# prepared cross-validation dataset
DATA_FILE <- here::here("data/data_prepared", paste0(SITE_NAME, "_cv.rds"))

# predictor set (columns absent from the data are silently dropped)
FEATURE_SET <- c(
  # --- Radiation ---
  "PPFD",                # photosynthetic photon flux density (µmol m-2 s-1)
  "Rg",                  # global shortwave radiation         (W m-2)
  # --- Atmosphere ---
  "VPD",                 # vapour pressure deficit            (hPa or kPa)
  "RH",                  # relative humidity                  (%)
  "Temp",                # air temperature                    (°C)
  "rain",                # precipitation per half-hour        (mm)
  "rain_rolling_24",     # cumulative rainfall, past 24 h     (mm)
  # --- Cyclical time encodings ---
  "hour_sin", "hour_cos", "doy_sin", "doy_cos", "month_sin", "month_cos",
  # --- Binary flag ---
  "night",
  # --- Management-specific (remove for unmanaged sites) ---
  "Grazing_days_since",  # activates the Phytomass Index (PI)
  "Fertiliser_days_since",
  "grass_height",
  "grass_biomass"
)

# TRUE = use the pre-computed PI_{gap_label} columns (Script 09) as predictors
PI_ENABLED <- TRUE

# output directory: results/{SITE}/{MANAGEMENT}/{MODEL}
RESULTS_DIR <- here::here("results", SITE_NAME, MANAGEMENT_CONDITION, MODEL_CHOICE)

# directory holding the model scripts
MODEL_SCRIPTS_DIR <- here::here("models")


# attach only the packages used across both scripts
suppressPackageStartupMessages({
  library(here); library(dplyr); library(glue); library(tibble); library(ranger)
})

# resolve the model script path and check both inputs exist
model_script <- file.path(MODEL_SCRIPTS_DIR, paste0(MODEL_CHOICE, "_CV.R"))
if (!file.exists(model_script)) stop("Model script not found: ", model_script)
if (!file.exists(DATA_FILE))    stop("Data file not found: ", DATA_FILE)

# load the prepared data
df <- readRDS(DATA_FILE)

# keep only the FEATURE_SET columns present in the data
predictors <- FEATURE_SET[FEATURE_SET %in% names(df)]
if (!length(predictors)) stop("None of the FEATURE_SET columns are in the data.")

# create the results directory (no-op if it already exists)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# pass settings to RF_CV.R via its own environment, then source it
run_env <- new.env(parent = globalenv())
run_env$df          <- df
run_env$predictors  <- predictors
run_env$RESULTS_DIR <- RESULTS_DIR
run_env$PI_ENABLED  <- PI_ENABLED

source(model_script, local = run_env)