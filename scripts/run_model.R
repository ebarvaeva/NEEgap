# =============================================================================
# run_model.R — Gap-Filling Entry Point: RF | MLP | XGBoost
# =============================================================================
#
# PURPOSE
#   This script is the single entry point for running any of the three
#   machine-learning gap-filling models (Random Forest, Multi-Layer
#   Perceptron, or XGBoost) on a half-hourly eddy-covariance flux dataset.
#   It gap-fills three carbon flux target variables simultaneously:
#
#     NEE  — Net Ecosystem Exchange        (µmol CO₂ m⁻² s⁻¹)
#     Reco — Ecosystem Respiration         (µmol CO₂ m⁻² s⁻¹)
#     GPP  — Gross Primary Production      (µmol CO₂ m⁻² s⁻¹)
#
#   Reco and GPP are not gap-filled directly; they are derived from the
#   gap-filled NEE series through flux partitioning (Lloyd-Taylor temperature
#   response for Reco; non-rectangular hyperbolic light response for GPP).
#   Partitioning parameters are estimated separately for each artificial gap
#   and each regrowth period, ensuring independence between training and
#   prediction.
#
# INPUTS
#   data/{SITE_NAME}.rds
#     A prepared half-hourly data frame (output of data_preparation/07_bind_years.R).
#     Required columns:
#       timestamp          — POSIXct, UTC, 30-min intervals
#       NEE_orig           — measured NEE (µmol m⁻² s⁻¹); NA where gapped
#       PPFD               — photosynthetic photon flux density
#       Temp               — air temperature (°C)
#       night              — binary night flag (1 = PPFD < 10)
#       + all predictor columns listed in FEATURE_SET
#
# OUTPUTS  (written to results/{SITE_NAME}/{MANAGEMENT_CONDITION}/{MODEL_CHOICE}/)
#   df_cv_all_predictions.rds
#     The original data frame augmented with 12 prediction columns:
#       {NEE|Reco|GPP}_{S|M|L|VL}_{model}_predicted
#
#   run_info.txt          — run metadata (timestamp, site, model, predictors)
#   progress.log          — timestamped progress messages
#
#   RF only:
#     variable_importance/
#       rf_variable_importance_{SIZE}_all.rds   — stacked per-gap VI table
#       rf_variable_importance_{SIZE}_all.csv
#       rf_vi_{SIZE}_{LABEL}.rds                — individual gap VI objects
#
#   MLP only:
#     training_loss_curves/{TARGET}/
#       MLP_{TARGET}_{SIZE}_{LABEL}_{loss|mae}.png
#
#   XGBoost only:
#     boosting_evaluation_curves/{TARGET}/
#       XGB_{TARGET}_{SIZE}_{LABEL}_rmse.png
#
# SECTIONS
#   1. User Settings     — site name, management condition, model, predictor set
#   2. Load Libraries    — model-specific packages loaded conditionally
#   3. Validate Settings — checks script and data file existence
#   4. Load Data         — reads the prepared RDS and intersects predictor set
#   5. Create Output Dir — creates results folder and logs run configuration
#   6. Source Model      — passes data to the selected model CV script
#
# HOW TO RUN
#   1. Open this script with the project root set as the working directory.
#   2. Edit the USER SETTINGS block (Section 1) for your site and model.
#   3. Source the script:  source("scripts/run_model.R")
#
# ADAPTING FOR A NEW SITE
#   1. Change SITE_NAME and MANAGEMENT_CONDITION.
#   2. Point DATA_FILE to the corresponding .rds file.
#   3. Set MODEL_CHOICE to "RF", "MLP", or "XGBoost".
#   4. Review FEATURE_SET and remove columns absent from the new site's data.
#
# =============================================================================


# =============================================================================
# SECTION 1 — User Settings
# =============================================================================
# Edit the variables below before running.  All other sections are automatic.

# --- Site and experimental condition -----------------------------------------
SITE_NAME            <- "JC1"       # Site identifier; used to locate data and name output directories
MANAGEMENT_CONDITION <- "managed"   # "managed" (with grazing/fertilisation data) or "unmanaged"

# --- Model choice ------------------------------------------------------------
MODEL_CHOICE <- "RF"   # "RF"  |  "MLP"  |  "XGBoost"

# --- Input data file ---------------------------------------------------------
DATA_FILE <- here::here(
  "data",
  paste0(SITE_NAME, ".rds")
)

# --- Predictor (feature) set -------------------------------------------------
# Columns listed here but absent from the data are silently dropped at runtime.
# Including "Grazing_days_since" activates two additional mechanisms:
#   (1) The Phytomass Index (PI) is computed per artificial gap (Section 5
#       of the model scripts) and appended to the feature set automatically.
#   (2) Flux partitioning for Reco and GPP is performed per regrowth period
#       (defined by grazing events) rather than globally.
# Remove management-specific columns for unmanaged sites.
FEATURE_SET <- c(
  # --- Radiation ---
  "PPFD",                    # photosynthetic photon flux density (µmol m⁻² s⁻¹)
  "Rg",                      # global shortwave radiation         (W m⁻²)
  # --- Atmosphere ---
  "VPD",                     # vapour pressure deficit            (hPa or kPa)
  "RH",                      # relative humidity                  (%)
  "Temp",                    # air temperature                    (°C)
  "rain",                    # precipitation per half-hour        (mm)
  "rain_rolling_24",         # cumulative rainfall, past 24 h     (mm)
  # --- Cyclical time encodings (sine/cosine pairs remove circular discontinuities) ---
  "hour_sin",  "hour_cos",   # hour of day
  "doy_sin",   "doy_cos",    # day of year
  "month_sin", "month_cos",  # month of year
  # --- Season dummies ---
  "Winter", "Spring", "Summer", "Autumn",
  # --- Binary flag ---
  "night",                   # 1 = PPFD < 10 µmol m⁻² s⁻¹
  # --- Management-specific predictors (remove for unmanaged sites) ---
  "Grazing_days_since",      # days since last grazing event; activates PI
  "Fertiliser_days_since",   # days since last fertiliser application
  "N",                       # cumulative N applied (kg N ha⁻¹)
  "grass_height",            # interpolated grass height (cm)
  "grass_biomass"            # interpolated dry biomass (kg DM ha⁻¹)
)

# --- Output directory --------------------------------------------------------
RESULTS_BASE_DIR <- here::here("results")
RESULTS_DIR      <- file.path(RESULTS_BASE_DIR,
                               SITE_NAME,
                               MANAGEMENT_CONDITION,
                               MODEL_CHOICE)

# --- Model script path -------------------------------------------------------
MODEL_SCRIPTS_DIR <- here::here("models")

# =============================================================================
# END OF USER SETTINGS — do not edit below this line for routine use
# =============================================================================


# =============================================================================
# SECTION 2 — Load Libraries
# =============================================================================
# Core packages are loaded unconditionally; model-specific packages are loaded
# only when required to avoid unnecessary dependency warnings.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(lubridate); library(rlang)
  library(glue); library(purrr); library(tidyr); library(ggplot2)
  library(tibble); library(zoo)
})

if (MODEL_CHOICE == "RF")      library(ranger)
if (MODEL_CHOICE == "XGBoost") library(xgboost)
if (MODEL_CHOICE == "MLP") {
  library(reticulate); library(tensorflow); library(keras3)
}


# =============================================================================
# SECTION 3 — Validate Settings
# =============================================================================
# Ensures the model script and data file exist before any computation begins,
# providing clear error messages if either is missing.

if (!MODEL_CHOICE %in% c("RF","MLP","XGBoost"))
  stop("MODEL_CHOICE must be 'RF', 'MLP', or 'XGBoost'.")

model_script <- file.path(MODEL_SCRIPTS_DIR, paste0(MODEL_CHOICE, "_CV.R"))
if (!file.exists(model_script))
  stop("Model script not found: ", model_script)
if (!file.exists(DATA_FILE))
  stop("Data file not found: ", DATA_FILE)


# =============================================================================
# SECTION 4 — Load Data and Build Predictor List
# =============================================================================
# The prepared RDS is loaded and the predictor set is intersected with the
# columns actually present in the data.  This allows FEATURE_SET to be
# specified generously; missing columns are reported but do not cause errors.

message("Loading data: ", DATA_FILE)
df         <- readRDS(DATA_FILE)
rds_name   <- basename(DATA_FILE)

# Retain only those predictors present in the data
predictors <- FEATURE_SET[FEATURE_SET %in% names(df)]
if (!length(predictors))
  stop("None of the FEATURE_SET columns found in the data frame.")

message("Predictors available (", length(predictors), " of ", length(FEATURE_SET), "):\n  ",
        paste(predictors, collapse = ", "))

absent <- setdiff(FEATURE_SET, names(df))
if (length(absent))
  message("NOTE — absent predictors (will be skipped):\n  ", paste(absent, collapse=", "))

# Warn if Reco/GPP partitioning will be limited
if (!all(c("PPFD","Temp") %in% names(df)))
  message("NOTE — 'PPFD' or 'Temp' not found in data: Reco/GPP gap-filling will be SKIPPED.")
if (!"Grazing_days_since" %in% names(df))
  message("NOTE — 'Grazing_days_since' absent: Reco/GPP will use ONE GLOBAL regrowth period.")


# =============================================================================
# SECTION 5 — Create Output Directory and Log Run Settings
# =============================================================================
# The output directory is created if it does not already exist.  Run metadata
# is printed to the console for reference.

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

message("\n--- Run configuration ---")
message("  Site               : ", SITE_NAME)
message("  Management         : ", MANAGEMENT_CONDITION)
message("  Model              : ", MODEL_CHOICE)
message("  Data rows          : ", nrow(df))
message("  Results directory  : ", RESULTS_DIR)
message("  Model script       : ", model_script)
message("-------------------------\n")


# =============================================================================
# SECTION 6 — Source the Model Script
# =============================================================================
# The selected model CV script (RF_CV.R, MLP_CV.R, or XGBoost_CV.R) is
# sourced in an isolated environment.  The following variables are passed:
#
#   df          — the prepared data frame
#   predictors  — the active predictor column names
#   RESULTS_DIR — the output directory path
#   rds_name    — the source file name (recorded in run_info.txt)
#
# Executing in a new environment prevents the model script from modifying
# the global workspace and allows multiple runs to be sourced in sequence
# without variable conflicts.

run_env <- new.env(parent = globalenv())
run_env$df          <- df
run_env$predictors  <- predictors
run_env$RESULTS_DIR <- RESULTS_DIR
run_env$rds_name    <- rds_name

run_status <- tryCatch({
  source(model_script, local = run_env)
  "completed successfully"
}, error = function(e) paste("ERROR —", conditionMessage(e)))

message("\n=== Run status: ", run_status, " ===")
message("Results in: ", RESULTS_DIR)

# =========================== end run_model.R ==================================
