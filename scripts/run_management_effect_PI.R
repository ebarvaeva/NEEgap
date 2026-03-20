# =============================================================================
# run_management_effect_PI.R — Phytomass Index (PI) Ablation Study
# =============================================================================
#
# PURPOSE
#   Runs all three gap-filling models (RF, MLP, XGBoost) using the BASE
#   meteorological predictor set augmented with the Phytomass Index (PI),
#   but without including Grazing_days_since as a direct split feature.
#
#   The Phytomass Index is a derived predictor that captures the state of
#   sward regrowth between grazing events.  It is computed internally by
#   each model script (via add_phytomass_index()) on a per-artificial-gap
#   basis: for each gap label, PI is estimated from the night-time and
#   daytime NEE observations outside that gap, reflecting the net carbon
#   balance of the canopy in the period surrounding the gap.
#
#   Because PI is computed fresh for each gap from observations outside it,
#   its construction is consistent with the cross-validation design and
#   does not leak future information into the training set.
#
# DESIGN RATIONALE
#   The Phytomass Index requires Grazing_days_since to be present in the
#   data (to identify regrowth periods), but it must NOT be passed as a
#   predictor — doing so would conflate the effect of PI with the direct
#   effect of Grazing_days_since.  The separation is achieved as follows:
#
#     predictors = BASE only       (no Grazing_days_since as a split feature)
#     PI_ENABLED = TRUE            (injected into the run environment)
#
#   The model scripts contain a guard:
#     if (!exists("PI_ENABLED")) {
#       PI_ENABLED <- "Grazing_days_since" %in% predictors
#     }
#   When PI_ENABLED is pre-set here to TRUE, the guard leaves it untouched.
#   When other run scripts do not pre-set it, the model derives it from the
#   predictor list as before, preserving backward compatibility.
#
# INPUTS
#   data/data_prepared/{SITE}.rds   (JC1 and JC2)
#     Must contain Grazing_days_since in the data (not as predictor).
#     Sites without it are skipped with a descriptive message.
#
# OUTPUTS  (one folder per site x model)
#   results/{SITE}/management_effect/{MODEL}/phytomass_index/
#     df_cv_all_predictions.rds
#     variable_importance/   (RF only)
#     run_info.txt, progress.log
#
#   Total runs: 2 sites x 3 models = 6 runs.
#
# SECTIONS
#   1. Configuration   — model scripts, predictor set, PI data column
#   2. Loop            — iterates over sites and models
#
# HOW TO RUN
#   source("scripts/run_management_effect_PI.R")
#
# =============================================================================


# =============================================================================
# SECTION 1 — Configuration
# =============================================================================

suppressPackageStartupMessages({
  library(here);      library(dplyr);    library(lubridate); library(rlang)
  library(glue);      library(purrr);    library(tidyr);     library(ggplot2)
  library(tibble);    library(zoo)
  # Model-specific libraries are loaded by the sourced model scripts
})

MODEL_SCRIPTS_DIR <- here::here("models")

# BASE predictor set — strictly meteorological and temporal variables.
# Grazing_days_since is intentionally absent from this list; PI is activated
# via PI_ENABLED below rather than via the predictor list.
BASE_PREDICTORS <- c(
  "PPFD", "Rg",
  "VPD", "RH", "Temp", "rain", "rain_rolling_24",
  "hour_sin", "hour_cos",
  "doy_sin",  "doy_cos",
  "month_sin","month_cos",
  "Winter", "Spring", "Summer", "Autumn",
  "night"
)

# Data column required for PI computation — not added to predictors.
# The model script's add_phytomass_index() uses this column to identify
# regrowth periods; the column must be present in the data even though
# it will not appear in the feature matrix passed to the model.
PI_DATA_COLUMN <- "Grazing_days_since"

# PI-specific model scripts — these contain the PI_ENABLED guard.
MODELS <- list(
  RF      = "RF_PI_CV.R",
  MLP     = "MLP_PI_CV.R",
  XGBoost = "XGBoost_PI_CV.R"
)

SITES         <- c("JC1", "JC2")
OUTPUT_SUBDIR <- "phytomass_index"   # subfolder inside management_effect/{MODEL}/


# =============================================================================
# SECTION 2 — Loop: site x model
# =============================================================================
# For each site, Grazing_days_since is verified to be present in the data
# (required by add_phytomass_index() even though it is not a predictor).
# The predictor set is restricted to BASE only.  PI_ENABLED = TRUE is
# injected so the model computes the phytomass index per gap.

for (site in SITES) {

  data_file <- here::here("data", "data_prepared", paste0(site, ".rds"))
  if (!file.exists(data_file)) {
    message("Data not found, skipping: ", data_file); next
  }

  df       <- readRDS(data_file)
  rds_name <- basename(data_file)

  # Guard: PI computation requires Grazing_days_since in the raw data
  if (!PI_DATA_COLUMN %in% names(df)) {
    message(
      "\n[SKIP] Site: ", site,
      " -- '", PI_DATA_COLUMN, "' column not found in data.",
      "\n       add_phytomass_index() needs this column to compute PI."
    )
    next
  }

  # Predictor set: BASE only; Grazing_days_since is excluded as a split feature
  predictors <- BASE_PREDICTORS[BASE_PREDICTORS %in% names(df)]

  message("\n", strrep("=", 60))
  message("Site      : ", site)
  message("Predictors: BASE only (", length(predictors), " vars)")
  message("PI_ENABLED: TRUE  [injected -- not derived from predictor list]")
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

    # Isolated environment — BASE predictors only; PI_ENABLED forced to TRUE.
    # The model script will compute PI_{gap_label} columns internally and
    # append them to the feature set for each gap training run.
    run_env             <- new.env(parent = globalenv())
    run_env$df          <- df
    run_env$predictors  <- predictors   # BASE only
    run_env$RESULTS_DIR <- RESULTS_DIR
    run_env$rds_name    <- rds_name
    run_env$PI_ENABLED  <- TRUE         # activates per-gap phytomass index

    run_status <- tryCatch({
      source(model_script, local = run_env)
      "completed successfully"
    }, error = function(e) paste("ERROR --", conditionMessage(e)))

    message("  Status: ", run_status)
  }
}

message("\n", strrep("=", 60))
message("Phytomass Index management effect runs complete.")
message(strrep("=", 60))
