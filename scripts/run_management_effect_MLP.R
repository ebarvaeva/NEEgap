# =============================================================================
# run_management_effect_MLP.R — Management Variable Ablation Study: MLP
# =============================================================================
#
# PURPOSE
#   Identical in design to run_management_effect_RF.R, but uses the
#   Multi-Layer Perceptron (MLP) neural network as the gap-filling model.
#   Each management variable is added individually to the BASE meteorological
#   predictor set and the model is trained and evaluated in full cross-
#   validation for both sites.
#
#   Running this script alongside its RF and XGBoost counterparts allows
#   the contribution of each management variable to be assessed across
#   three structurally different model families, providing evidence of
#   generalisability rather than model-specific artefacts.
#
# INPUTS
#   data/data_prepared/{SITE}.rds   (JC1 and JC2)
#
# OUTPUTS  (one folder per site × management variable)
#   results/{SITE}/management_effect/MLP/{MGMT_VAR}/
#     df_cv_all_predictions.rds
#     training_loss_curves/{TARGET}/   — per-gap training and validation loss plots
#     run_info.txt, progress.log
#
#   Total runs: 2 sites × 5 management variables = 10 runs.
#
# BASE PREDICTOR SET
#   Meteorological:  PPFD, Rg, VPD, RH, Temp, rain, rain_rolling_24
#   Temporal:        hour_sin/cos, doy_sin/cos, month_sin/cos
#   Season dummies:  Winter, Spring, Summer, Autumn
#   Binary flag:     night
#
# MANAGEMENT VARIABLES TESTED (one at a time)
#   Grazing_days_since    — days since last grazing event
#   Fertiliser_days_since — days since last fertilisation event
#   N                     — cumulative nitrogen applied (kg N ha⁻¹)
#   grass_height          — interpolated sward height (cm)
#   grass_biomass         — interpolated dry biomass (kg DM ha⁻¹)
#
# SECTIONS
#   1. Configuration   — model path, predictor sets, sites
#   2. Loop            — iterates over sites and management variables
#
# REQUIREMENTS
#   Python environment with TensorFlow/Keras accessible via reticulate.
#   The Python environment must be configured before running.
#
# HOW TO RUN
#   source("scripts/run_management_effect_MLP.R")
#
# =============================================================================


# =============================================================================
# SECTION 1 — Configuration
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(lubridate); library(rlang)
  library(glue); library(purrr); library(tidyr); library(ggplot2)
  library(tibble); library(zoo); library(reticulate); library(tensorflow); library(keras3)
})

MODEL_CHOICE      <- "MLP"
MODEL_SCRIPTS_DIR <- here::here("models")
model_script      <- file.path(MODEL_SCRIPTS_DIR, "MLP_CV.R")
if (!file.exists(model_script)) stop("Model script not found: ", model_script)

# BASE predictor set — no management variables included.
BASE_PREDICTORS <- c(
  "PPFD", "Rg",
  "VPD", "RH", "Temp", "rain", "rain_rolling_24",
  "hour_sin", "hour_cos",
  "doy_sin",  "doy_cos",
  "month_sin","month_cos",
  "Winter", "Spring", "Summer", "Autumn",
  "night"
)

# Management variables to test one at a time.
MGMT_VARS <- c(
  "Grazing_days_since",
  "Fertiliser_days_since",
  "N",
  "grass_height",
  "grass_biomass"
)

SITES <- c("JC1", "JC2")


# =============================================================================
# SECTION 2 — Loop: site × management variable
# =============================================================================
# For each site, the dataset is loaded once.  The inner loop adds one
# management variable at a time, creates the output directory, and sources
# MLP_CV.R in an isolated environment.

for (site in SITES) {

  data_file <- here::here("data", "data_prepared", paste0(site, ".rds"))
  if (!file.exists(data_file)) {
    message("Data not found, skipping: ", data_file); next
  }
  df       <- readRDS(data_file)
  rds_name <- basename(data_file)
  message("\n", strrep("=", 60))
  message("Site: ", site, "  |  Model: ", MODEL_CHOICE)
  message(strrep("=", 60))

  for (mgmt_var in MGMT_VARS) {

    if (!mgmt_var %in% names(df)) {
      message("  Skipping ", mgmt_var, " — not found in ", site, " data.")
      next
    }

    feature_set <- c(BASE_PREDICTORS, mgmt_var)
    predictors  <- feature_set[feature_set %in% names(df)]

    RESULTS_DIR <- here::here("results", site, "management_effect",
                               MODEL_CHOICE, mgmt_var)
    dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

    message("\n", strrep("-", 50))
    message("  Predictor set: BASE + ", mgmt_var)
    message("  Results dir  : ", RESULTS_DIR)
    message(strrep("-", 50))

    run_env             <- new.env(parent = globalenv())
    run_env$df          <- df
    run_env$predictors  <- predictors
    run_env$RESULTS_DIR <- RESULTS_DIR
    run_env$rds_name    <- rds_name

    run_status <- tryCatch({
      source(model_script, local = run_env)
      "completed successfully"
    }, error = function(e) paste("ERROR —", conditionMessage(e)))

    message("  Status: ", run_status)
  }
}

message("\n", strrep("=", 60))
message("Management effect MLP runs complete.")
message(strrep("=", 60))
