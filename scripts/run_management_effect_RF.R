# =============================================================================
# run_management_effect_RF.R — Management Variable Ablation Study: RF
# =============================================================================
#
# PURPOSE
#   Executes a systematic ablation study to quantify how each individual
#   management variable contributes to gap-filling performance when added
#   to the BASE meteorological predictor set.  The Random Forest model is
#   trained separately for each combination of site × management variable,
#   producing an isolated estimate of each variable's predictive value.
#
#   The BASE predictor set contains only meteorological and temporal
#   variables (no management information).  Each run adds exactly one
#   management variable to this base, allowing the improvement in gap-
#   filling accuracy to be attributed specifically to that variable.
#
#   This design mirrors a standard feature-addition experiment and produces
#   results directly comparable across variables and models.
#
# INPUTS
#   data/data_prepared/{SITE}.rds
#     Prepared dataset for JC1 and JC2 (output of data_preparation pipeline).
#     Must contain the BASE predictor columns and each management variable.
#
# OUTPUTS  (one folder per site × management variable combination)
#   results/{SITE}/management_effect/RF/{MGMT_VAR}/
#     df_cv_all_predictions.rds       — predictions for all targets × gap sizes
#     variable_importance/            — RF impurity importance per gap label
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
#   Grazing_days_since    — days elapsed since the last grazing event
#   Fertiliser_days_since — days elapsed since the last fertilisation event
#   N                     — cumulative nitrogen applied (kg N ha⁻¹)
#   grass_height          — interpolated sward height (cm)
#   grass_biomass         — interpolated dry biomass (kg DM ha⁻¹)
#
# NOTE ON GRAZING_DAYS_SINCE
#   When Grazing_days_since is included in the predictor set, the model
#   script automatically activates the Phytomass Index (PI) and per-regrowth
#   flux partitioning for Reco/GPP.  To isolate the effect of Grazing_days_since
#   without PI, use run_management_effect_GrazingDaysSince.R instead.
#
# SECTIONS
#   1. Configuration   — model path, predictor sets, sites
#   2. Loop            — iterates over sites and management variables
#
# HOW TO RUN
#   source("scripts/run_management_effect_RF.R")
#
# =============================================================================


# =============================================================================
# SECTION 1 — Configuration
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(lubridate); library(rlang)
  library(glue); library(purrr); library(tidyr); library(ggplot2)
  library(tibble); library(zoo); library(ranger)
})

MODEL_CHOICE      <- "RF"
MODEL_SCRIPTS_DIR <- here::here("models")
model_script      <- file.path(MODEL_SCRIPTS_DIR, "RF_CV.R")
if (!file.exists(model_script)) stop("Model script not found: ", model_script)

# BASE predictor set — meteorological and temporal variables only.
# No management variables are included here; each is added individually in the loop.
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
# Each variable is appended to BASE_PREDICTORS for one isolated run.
MGMT_VARS <- c(
  "Grazing_days_since",
  "Fertiliser_days_since",
  "N",
  "grass_height",
  "grass_biomass"
)

SITES <- c("JC1", "JC2")   # Both sites are processed; change to run one at a time


# =============================================================================
# SECTION 2 — Loop: site × management variable
# =============================================================================
# For each site, the prepared dataset is loaded once.  The inner loop then
# constructs the predictor set for each management variable, creates the
# output directory, and sources the RF model script in an isolated environment.
# Variables absent from the site's data are skipped with a message.

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

    # Skip if this management variable is absent from the site's data
    if (!mgmt_var %in% names(df)) {
      message("  Skipping ", mgmt_var, " — not found in ", site, " data.")
      next
    }

    # Construct predictor set: BASE + this management variable only
    feature_set <- c(BASE_PREDICTORS, mgmt_var)
    predictors  <- feature_set[feature_set %in% names(df)]

    RESULTS_DIR <- here::here("results", site, "management_effect",
                               MODEL_CHOICE, mgmt_var)
    dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

    message("\n", strrep("-", 50))
    message("  Predictor set: BASE + ", mgmt_var)
    message("  Results dir  : ", RESULTS_DIR)
    message(strrep("-", 50))

    # Source in an isolated environment so that each run is fully independent
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
message("Management effect RF runs complete.")
message(strrep("=", 60))
