# =============================================================================
# RUN SCRIPT — Eddy Covariance Gap-Filling 
# =============================================================================
#
# PURPOSE
#   Entry point for running any of the three gap-filling models (RF, MLP,
#   XGBoost) on a flux-tower dataset.  This version gap-fills three target
#   variables: NEE, Reco, and GPP.
#
#   Reco and GPP are derived by flux partitioning (see models/RF_CV.R §7
#   for the full algorithm), using the Lloyd & Taylor temperature response
#   for Reco and a non-rectangular hyperbolic light response for GPP.
#   Partitioning is performed PER REGROWTH PERIOD and PER ARTIFICIAL GAP,
#   ensuring that no observations from inside a gap influence the parameter
#   estimates used to fill that same gap.
#
# REQUIRED INPUT DATA COLUMNS
#   timestamp        — POSIXct or parseable date-time string (half-hourly)
#   NEE_orig         — net ecosystem exchange                (µmol m⁻² s⁻¹)
#   PPFD             — photon flux density — REQUIRED for Reco/GPP
#   Temp             — air temperature     — REQUIRED for Reco/GPP
#   night            — binary flag (1 = night, 0 = day)
#   Grazing_days_since — days since last grazing  (optional; enables PI and
#                        per-regrowth partitioning; if absent, one global
#                        regrowth period is used)
#   + all additional columns listed in FEATURE_SET
#
# PHYTOMASS INDEX (PI)
#   When PI_ENABLED = TRUE, each model receives PI_{gap_label} as an
#   additional predictor (e.g. PI_S1 for gap S1, PI_M3 for gap M3).
#   These columns must already exist in the input RDS, produced by
#   Script 09 (phytomass_index.R) — i.e. DATA_FILE should point to
#   JCi_cv.rds.  PI is computed exclusively from
#   outside-gap NEE_orig observations (Lohila et al. 2004), so it
#   carries no information from the test set into the model.
#   Set PI_ENABLED = FALSE to run without PI for any site or model.
#
# OUTPUTS PER MODEL
#   results/{SITE}/{MANAGEMENT}/{MODEL}/
#     df_cv_all_predictions.rds           ← 12 prediction columns (NEE/Reco/GPP × S/M/L/VL)
#     run_info.txt
#     progress.log
#
#     RF only:
#       variable_importance/NEE/          ← rf_variable_importance_{SIZE}_all.rds / .csv
#       variable_importance/Reco/
#       variable_importance/GPP/
#
#     MLP only:
#       training_loss_curves/NEE/         ← MLP_NEE_{SIZE}_{LABEL}_{loss|mae}.png
#       training_loss_curves/Reco/
#       training_loss_curves/GPP/
#
#     XGBoost only:
#       boosting_evaluation_curves/NEE/   ← XGB_NEE_{SIZE}_{LABEL}_rmse.png
#       boosting_evaluation_curves/Reco/
#       boosting_evaluation_curves/GPP/
#
# ADAPTING FOR A DIFFERENT SITE
#   1. Change SITE_NAME and MANAGEMENT_CONDITION.
#   2. Point DATA_FILE to the new .rds file (JCi_cv.rds if PI_ENABLED = TRUE).
#   3. Set MODEL_CHOICE to "RF", "MLP", or "XGBoost".
#   4. Set PI_ENABLED = TRUE / FALSE as needed.
#   5. Adjust FEATURE_SET if needed.
#
# =============================================================================




# =============================================================================
# USER SETTINGS — edit this block for each new site / run
# =============================================================================

# --- Site and experimental condition -----------------------------------------
SITE_NAME            <- "JC1"       # change to "JC2", "JC3", etc.
MANAGEMENT_CONDITION <- "managed"   # "managed"   or   "unmanaged"

# --- Model choice ------------------------------------------------------------
MODEL_CHOICE <- "RF"   # "RF"  |  "MLP"  |  "XGBoost"

# --- Input data file ---------------------------------------------------------
DATA_FILE <- here::here(
  "data",
  paste0(SITE_NAME, ".rds")
)

# --- Predictor (feature) set -------------------------------------------------
# Columns absent from the data are silently dropped.
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
  # --- Cyclical time encodings ---
  "hour_sin",  "hour_cos",   # hour of day
  "doy_sin",   "doy_cos",    # day of year
  "month_sin", "month_cos",  # month of year
  # --- Season dummies ---
  "Winter", "Spring", "Summer", "Autumn",
  # --- Binary flag ---
  "night",
  # --- Management-specific (remove for unmanaged sites) ---
  "Grazing_days_since",      # activates PI + per-regrowth partitioning
  "Fertiliser_days_since",
  "N",
  "grass_height",
  "grass_biomass"
)

# --- Phytomass Index ----------------------------------------------------------
# Set TRUE to use pre-computed PI_{gap_label} columns (from Script 09) as
# predictors.  Requires JCi_cv.rds to contain PI columns.
PI_ENABLED <- TRUE

# --- Output directory --------------------------------------------------------
RESULTS_BASE_DIR <- here::here("results")
RESULTS_DIR      <- file.path(RESULTS_BASE_DIR,
                               SITE_NAME,
                               MANAGEMENT_CONDITION,
                               MODEL_CHOICE)

# --- Model script path -------------------------------------------------------
MODEL_SCRIPTS_DIR <- here::here("models")

# =============================================================================
# END OF USER SETTINGS
# =============================================================================


# =============================================================================
# LOAD LIBRARIES
# =============================================================================

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
# VALIDATE SETTINGS
# =============================================================================

if (!MODEL_CHOICE %in% c("RF","MLP","XGBoost"))
  stop("MODEL_CHOICE must be 'RF', 'MLP', or 'XGBoost'.")

model_script <- file.path(MODEL_SCRIPTS_DIR, paste0(MODEL_CHOICE, "_CV.R"))
if (!file.exists(model_script))
  stop("Model script not found: ", model_script)
if (!file.exists(DATA_FILE))
  stop("Data file not found: ", DATA_FILE)


# =============================================================================
# LOAD DATA AND BUILD PREDICTOR LIST
# =============================================================================

message("Loading data: ", DATA_FILE)
df         <- readRDS(DATA_FILE)
rds_name   <- basename(DATA_FILE)

predictors <- FEATURE_SET[FEATURE_SET %in% names(df)]
if (!length(predictors))
  stop("None of the FEATURE_SET columns found in the data frame.")

message("Predictors available (", length(predictors), " of ", length(FEATURE_SET), "):\n  ",
        paste(predictors, collapse = ", "))

absent <- setdiff(FEATURE_SET, names(df))
if (length(absent))
  message("NOTE — absent predictors (will be skipped):\n  ", paste(absent, collapse=", "))

# Check Reco/GPP requirements
if (!all(c("PPFD","Temp") %in% names(df)))
  message("NOTE — 'PPFD' or 'Temp' not found in data: Reco/GPP gap-filling will be SKIPPED.")
if (!"Grazing_days_since" %in% names(df))
  message("NOTE — 'Grazing_days_since' absent: Reco/GPP will use ONE GLOBAL regrowth period.")


# =============================================================================
# CREATE OUTPUT DIRECTORY AND LOG RUN SETTINGS
# =============================================================================

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
# SOURCE THE MODEL SCRIPT
# =============================================================================

run_env <- new.env(parent = globalenv())
run_env$df          <- df
run_env$predictors  <- predictors
run_env$RESULTS_DIR <- RESULTS_DIR
run_env$rds_name    <- rds_name
run_env$PI_ENABLED  <- PI_ENABLED 

run_status <- tryCatch({
  source(model_script, local = run_env)
  "completed successfully"
}, error = function(e) paste("ERROR —", conditionMessage(e)))

message("\n=== Run status: ", run_status, " ===")
message("Results in: ", RESULTS_DIR)

# =========================== end run_model.R (v2) =============================
