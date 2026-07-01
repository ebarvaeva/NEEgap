# run_management_effect_PI.R — Phytomass Index batch runner (all models)
# Runs RF, XGBoost, and MLP once per site with the BASE predictor set and
# PI_ENABLED = TRUE, so each model script appends the pre-computed PI_{gap_label}
# column (Script 09) for the gap it is filling. Compare against the BASE-only
# runs to isolate the contribution of the Phytomass Index. RF uses no season
# dummies (seasonality is captured by the sin/cos encodings; hard 0/1 boundaries
# hurt trees); XGBoost and MLP include them.
# Required (in each {SITE}_cv.rds, from Scripts 08-09):
#   timestamp, NEE_orig, night, gap-flag + masked-NEE columns, BASE predictors,
#   and PI_{gap_label} columns from Script 09 (the site is skipped if absent)
# Writes: nothing itself — each models/{MODEL}_CV.R saves predictions to
#   results/{SITE}/management_effect/{MODEL}/phytomass_index/
# To adapt: edit SITES and MODELS below (currently SITES = JC2; MODELS run in
#   the order RF, XGBoost, MLP).

# Packages (model-specific libraries are loaded on first use in the loop)
suppressPackageStartupMessages({
  library(here); library(dplyr); library(lubridate); library(rlang)
  library(glue); library(purrr); library(tidyr); library(ggplot2)
  library(tibble); library(zoo)
})

MODEL_SCRIPTS_DIR <- here::here("models")

# RF BASE set — no season dummies
BASE_PREDICTORS_RF <- c(
  "PPFD", "Rg",
  "VPD", "RH", "Temp", "rain", "rain_rolling_24",
  "hour_sin", "hour_cos",
  "doy_sin",  "doy_cos",
  "month_sin","month_cos",
  "night"
)

# XGBoost & MLP BASE set — RF set plus season dummies
BASE_PREDICTORS_LINEAR <- c(
  BASE_PREDICTORS_RF,
  "Winter", "Spring", "Summer", "Autumn"
)

SITES  <- c("JC1", "JC2")                    # sites to run 
MODELS <- c("RF", "XGBoost", "MLP")   # run in this order

# For each site: load data, verify PI columns, then run every model with PI
for (site in SITES) {

  # {site}_cv.rds holds gap-flag, masked-NEE, and PI columns (Scripts 08-09)
  data_file <- here::here("data", "data_prepared", paste0(site, "_cv.rds"))
  if (!file.exists(data_file)) next
  df       <- readRDS(data_file)
  rds_name <- basename(data_file)

  # Skip this site if no PI_ columns are present (Script 09 must run first)
  pi_cols <- names(df)[grepl("^PI_", names(df))]
  if (!length(pi_cols)) next

  for (model_choice in MODELS) {

    # RF drops season dummies; XGBoost and MLP keep them
    base_set   <- if (model_choice == "RF") BASE_PREDICTORS_RF else BASE_PREDICTORS_LINEAR
    predictors <- base_set[base_set %in% names(df)]

    # Load model-specific libraries on first use
    if (model_choice == "RF")      library(ranger)
    if (model_choice == "XGBoost") library(xgboost)
    if (model_choice == "MLP") {
      library(reticulate); library(tensorflow); library(keras3)
    }

    # Locate this model's CV script; skip if missing
    model_script <- file.path(MODEL_SCRIPTS_DIR, paste0(model_choice, "_CV.R"))
    if (!file.exists(model_script)) next

    # Output folder for this model's PI run
    RESULTS_DIR <- here::here("results", site, "management_effect", model_choice, "phytomass_index")
    dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

    # Pass data + settings to the model script in an isolated environment
    run_env             <- new.env(parent = globalenv())
    run_env$df          <- df
    run_env$predictors  <- predictors
    run_env$RESULTS_DIR <- RESULTS_DIR
    run_env$rds_name    <- rds_name
    run_env$PI_ENABLED  <- TRUE   # PI_{gap_label} appended per gap inside model script

    # Run the model (an error in one run does not stop the batch)
    tryCatch(source(model_script, local = run_env), error = function(e) NULL)
  }
}
