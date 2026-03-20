# =============================================================================
# run_management_effect_GrazingDaysSince.R — Grazing Days Since (No PI)
# =============================================================================
#
# PURPOSE
#   Runs all three gap-filling models (RF, MLP, XGBoost) with the predictor
#   set BASE + Grazing_days_since, where Grazing_days_since is treated as
#   a direct numeric split feature only — without activating the Phytomass
#   Index (PI).
#
#   This run is distinct from the standard management effect scripts
#   (run_management_effect_RF.R etc.) because those scripts automatically
#   activate PI whenever Grazing_days_since appears in the predictor list
#   (PI_ENABLED is derived from "Grazing_days_since" %in% predictors in each
#   model script).  Here, PI_ENABLED is explicitly set to FALSE in each
#   run environment, so the model sees Grazing_days_since as a raw predictor
#   but does not compute the per-gap phytomass index.
#
#   This isolation allows the contribution of Grazing_days_since alone to be
#   separated from the contribution of the PI that it indirectly enables,
#   which is quantified separately by run_management_effect_PI.R.
#
# INPUTS
#   data/data_prepared/{SITE}.rds   (JC1 and JC2)
#     Must contain Grazing_days_since; sites without it are skipped.
#
# OUTPUTS  (one folder per site × model)
#   results/{SITE}/management_effect/{MODEL}/Grazing_days_since/
#     df_cv_all_predictions.rds
#     variable_importance/   (RF only)
#     run_info.txt, progress.log
#
#   Total runs: 2 sites × 3 models = 6 runs.
#
# SECTIONS
#   1. Configuration   — model scripts, predictor set, sites
#   2. Loop            — iterates over sites and models
#
# HOW TO RUN
#   source("scripts/run_management_effect_GrazingDaysSince.R")
#
# =============================================================================


# =============================================================================
# SECTION 1 — Configuration
# =============================================================================

suppressPackageStartupMessages({
  library(here);      library(dplyr);    library(lubridate); library(rlang)
  library(glue);      library(purrr);    library(tidyr);     library(ggplot2)
  library(tibble);    library(zoo)
})

MODEL_SCRIPTS_DIR <- here::here("models")

# BASE predictor set — meteorological and temporal variables only.
BASE_PREDICTORS <- c(
  "PPFD", "Rg",
  "VPD", "RH", "Temp", "rain", "rain_rolling_24",
  "hour_sin", "hour_cos",
  "doy_sin",  "doy_cos",
  "month_sin","month_cos",
  "Winter", "Spring", "Summer", "Autumn",
  "night"
)

# The management variable under investigation in this script.
# It is added to BASE predictors as a direct numeric feature.
MGMT_VAR <- "Grazing_days_since"

# Model scripts: each model has a dedicated variant without PI logic.
# These scripts do not contain the PI activation guard, ensuring that
# PI_ENABLED = FALSE passed via run_env is respected unconditionally.
MODELS <- list(
  RF      = "RF_Grazing_CV.R",
  MLP     = "MLP_Grazing_CV.R",
  XGBoost = "XGBoost_Grazing_CV.R"
)

SITES         <- c("JC1", "JC2")
OUTPUT_SUBDIR <- "Grazing_days_since"   # subfolder name inside management_effect/{MODEL}/


# =============================================================================
# SECTION 2 — Loop: site × model
# =============================================================================
# For each site, the dataset is loaded once and Grazing_days_since is
# verified to be present.  The predictor set is then constructed as
# BASE + Grazing_days_since.  PI_ENABLED = FALSE is injected explicitly
# into each run environment so that the model script does not compute
# the phytomass index even though Grazing_days_since is in predictors.

for (site in SITES) {

  data_file <- here::here("data", "data_prepared", paste0(site, ".rds"))
  if (!file.exists(data_file)) {
    message("Data not found, skipping: ", data_file); next
  }

  df       <- readRDS(data_file)
  rds_name <- basename(data_file)

  # Skip sites where Grazing_days_since is not available
  if (!MGMT_VAR %in% names(df)) {
    message("\n[SKIP] Site: ", site, " — '", MGMT_VAR, "' not found in data.")
    next
  }

  # Predictor set: BASE + Grazing_days_since; filtered to columns in this site's data
  predictors <- c(BASE_PREDICTORS, MGMT_VAR)
  predictors <- predictors[predictors %in% names(df)]

  message("\n", strrep("=", 60))
  message("Site      : ", site)
  message("Predictors: BASE + ", MGMT_VAR, " (", length(predictors), " vars)")
  message("PI_ENABLED: FALSE  [injected — PI suppressed despite Grazing_days_since in predictors]")
  message(strrep("=", 60))

  for (model_name in names(MODELS)) {

    model_script <- file.path(MODEL_SCRIPTS_DIR, MODELS[[model_name]])
    if (!file.exists(model_script)) {
      message("  [SKIP] Script not found: ", model_script); next
    }

    RESULTS_DIR <- here::here(
      "results", site, "management_effect", model_name, OUTPUT_SUBDIR
    )
    dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

    message("\n", strrep("-", 50))
    message("  Model      : ", model_name)
    message("  Results dir: ", RESULTS_DIR)
    message(strrep("-", 50))

    run_env             <- new.env(parent = globalenv())
    run_env$df          <- df
    run_env$predictors  <- predictors
    run_env$RESULTS_DIR <- RESULTS_DIR
    run_env$rds_name    <- rds_name
    run_env$PI_ENABLED  <- FALSE   # PI explicitly disabled; Grazing_days_since is a direct feature

    run_status <- tryCatch({
      source(model_script, local = run_env)
      "completed successfully"
    }, error = function(e) paste("ERROR —", conditionMessage(e)))

    message("  Status: ", run_status)
  }
}

message("\n", strrep("=", 60))
message("Grazing_days_since management effect runs complete.")
message(strrep("=", 60))
