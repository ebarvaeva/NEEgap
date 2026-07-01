# run_management_effect_XGBoost.R — XGBoost feature-addition batch runner
# Trains XGBoost once per management variable, each time on the BASE predictor
# set plus ONE management variable, to isolate that variable's marginal effect
# on NEE gap-filling. PI_ENABLED is always FALSE here (raw predictor effect
# only; use run_management_effect_PI.R for the Phytomass Index). Loops over
# every site in SITES and every variable in MGMT_VARS, sourcing
# models/XGBoost_CV.R in an isolated environment for each combination.
# Required (in each {SITE}_cv.rds, from Scripts 08-09):
#   timestamp — POSIXct / parseable half-hourly date-time
#   NEE_orig  — net ecosystem exchange (µmol m-2 s-1)
#   night     — binary flag (1 = night, 0 = day)
#   gap-flag + masked-NEE columns (Script 08)
#   + BASE predictors and each MGMT_VARS column present (absent ones skipped)
# Writes: nothing itself — models/XGBoost_CV.R saves predictions to
#   results/{SITE}/management_effect/XGBoost/{mgmt_var}/
# To adapt: edit SITES and MGMT_VARS below (currently SITES = JC2, 4 variables).

# Packages (core + XGBoost-specific)
suppressPackageStartupMessages({
  library(here); library(dplyr); library(lubridate); library(rlang)
  library(glue); library(purrr); library(tidyr); library(ggplot2)
  library(tibble); library(zoo); library(xgboost)
})

# Locate the XGBoost CV script
MODEL_CHOICE      <- "XGBoost"
MODEL_SCRIPTS_DIR <- here::here("models")
model_script      <- file.path(MODEL_SCRIPTS_DIR, "XGBoost_CV.R")
if (!file.exists(model_script)) stop("Model script not found: ", model_script)

# BASE predictor set — no management variables (season dummies included)
BASE_PREDICTORS <- c(
  "PPFD", "Rg",
  "VPD", "RH", "Temp", "rain", "rain_rolling_24",
  "hour_sin", "hour_cos",
  "doy_sin",  "doy_cos",
  "month_sin","month_cos",
  "Winter", "Spring", "Summer", "Autumn",
  "night"
)

# Management variables, tested one at a time
MGMT_VARS <- c(
  "Grazing_days_since",
  "Fertiliser_days_since",
  "grass_height",
  "grass_biomass"
)

# Sites to run (edit to add "JC1")
SITES <- c("JC1", "JC2")

# For each site: load its prepared data, then loop over management variables
for (site in SITES) {

  # {site}_cv.rds holds the gap-flag and masked-NEE columns (Scripts 08-09)
  data_file <- here::here("data", "data_prepared", paste0(site, "_cv.rds"))
  if (!file.exists(data_file)) next
  df       <- readRDS(data_file)
  rds_name <- basename(data_file)

  for (mgmt_var in MGMT_VARS) {
    # Skip a management variable that is absent from this site's data
    if (!mgmt_var %in% names(df)) next

    # BASE + this one management variable, keeping only columns that exist
    feature_set <- c(BASE_PREDICTORS, mgmt_var)
    predictors  <- feature_set[feature_set %in% names(df)]

    # Output folder for this variable's run
    RESULTS_DIR <- here::here("results", site, "management_effect",
                              MODEL_CHOICE, mgmt_var)
    dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

    # Pass data + settings to the model script in an isolated environment
    run_env             <- new.env(parent = globalenv())
    run_env$df          <- df
    run_env$predictors  <- predictors
    run_env$RESULTS_DIR <- RESULTS_DIR
    run_env$rds_name    <- rds_name
    run_env$PI_ENABLED  <- FALSE   # raw predictor effect only — no PI

    # Run the model (an error in one run does not stop the batch)
    tryCatch(source(model_script, local = run_env), error = function(e) NULL)
  }
}
