# =============================================================================
# RUN SCRIPT — Management Effect Experiment
# =============================================================================
#
# PURPOSE
#   Entry point for running a single management effect experiment.
#   Each run tests one management variable (or PI) added on top of the BASE
#   predictor set, allowing you to isolate the contribution of each variable.
#
#   Run this script once per management event per model per site.
#   Six management events are available:
#
#     1. "Grazing_days_since"    — days since last grazing event
#     2. "Fertiliser_days_since" — days since last fertiliser application
#     3. "N"                     — nitrogen application amount
#     4. "grass_height"          — sward height at measurement
#     5. "grass_biomass"         — above-ground biomass
#     6. "PI"                    — Phytomass Index (no extra predictor column;
#                                  PI_{gap_label} columns appended per gap
#                                  by the model script via PI_ENABLED = TRUE)
#
#   Events 1–5 add the named column to the BASE predictor set with
#   PI_ENABLED = FALSE (testing the raw predictor in isolation).
#   Event 6 (PI) uses the BASE predictor set only with PI_ENABLED = TRUE
#   (testing the pre-computed PI signal in isolation).
#
# PRE-REQUISITES
#   - Script 08 and Script 09 must have been run for the chosen site so that
#     the _cv.rds file exists and contains gap flag, masked-NEE, and PI columns.
#   - For event 6 (PI), Script 09 must have been run (PI columns required).
#
# OUTPUTS
#   results/{SITE}/management_effect/{MODEL}/{MGMT_EVENT}/
#     df_cv_all_predictions.rds
#     run_info.txt
#     progress.log
#     variable_importance/   (RF only)
#     boosting_evaluation_curves/  (XGBoost only)
#     training_loss_curves/        (MLP only)
#
# ADAPTING FOR A DIFFERENT RUN
#   1. Set SITE_NAME.
#   2. Set MODEL_CHOICE.
#   3. Set MGMT_EVENT to one of the six options listed above.
#   4. Source the script.
#
# =============================================================================


# =============================================================================
# USER SETTINGS — edit this block for each new run
# =============================================================================

# --- Site --------------------------------------------------------------------
SITE_NAME <- "JC1"       # "JC1" or "JC2"

# --- Model -------------------------------------------------------------------
MODEL_CHOICE <- "RF"     # "RF"  |  "XGBoost"  |  "MLP"

# --- Management event --------------------------------------------------------
# Choose ONE of the six options below (copy exactly as shown):
#
#   "Grazing_days_since"      event 1 — BASE + grazing predictor,  PI = FALSE
#   "Fertiliser_days_since"   event 2 — BASE + fertiliser predictor, PI = FALSE
#   "N"                       event 3 — BASE + N predictor,          PI = FALSE
#   "grass_height"            event 4 — BASE + height predictor,     PI = FALSE
#   "grass_biomass"           event 5 — BASE + biomass predictor,    PI = FALSE
#   "PI"                      event 6 — BASE only, PI_ENABLED = TRUE
#
MGMT_EVENT <- "Grazing_days_since"

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

valid_events <- c("Grazing_days_since", "Fertiliser_days_since",
                  "N", "grass_height", "grass_biomass", "PI")
if (!MGMT_EVENT %in% valid_events)
  stop("MGMT_EVENT must be one of: ", paste(valid_events, collapse = ", "))

if (!MODEL_CHOICE %in% c("RF","MLP","XGBoost"))
  stop("MODEL_CHOICE must be 'RF', 'MLP', or 'XGBoost'.")

MODEL_SCRIPTS_DIR <- here::here("models")
model_script <- file.path(MODEL_SCRIPTS_DIR, paste0(MODEL_CHOICE, "_CV.R"))
if (!file.exists(model_script))
  stop("Model script not found: ", model_script)

data_file <- here::here("data", "data_prepared", paste0(SITE_NAME, "_cv.rds"))
if (!file.exists(data_file))
  stop("Data file not found: ", data_file,
       "\n  Run Scripts 08 and 09 first to generate the _cv.rds file.")


# =============================================================================
# BASE PREDICTOR SET
# =============================================================================
# Management columns are NOT included here — they are added one at a time
# below based on MGMT_EVENT.  PI columns are never listed here; they are
# appended automatically per gap inside the model script when PI_ENABLED = TRUE.

BASE_PREDICTORS <- c(
  # --- Radiation ---
  "PPFD", "Rg",
  # --- Atmosphere ---
  "VPD", "RH", "Temp", "rain", "rain_rolling_24",
  # --- Cyclical time encodings ---
  "hour_sin", "hour_cos",
  "doy_sin",  "doy_cos",
  "month_sin","month_cos",
  # --- Season dummies ---
  "Winter", "Spring", "Summer", "Autumn",
  # --- Binary flag ---
  "night"
)


# =============================================================================
# RESOLVE PREDICTOR SET AND PI FLAG FOR CHOSEN EVENT
# =============================================================================

if (MGMT_EVENT == "PI") {
  # Event 6: BASE only + PI signal injected per gap by the model script
  feature_set <- BASE_PREDICTORS
  PI_ENABLED  <- TRUE
} else {
  # Events 1–5: BASE + one management column, no PI
  feature_set <- c(BASE_PREDICTORS, MGMT_EVENT)
  PI_ENABLED  <- FALSE
}


# =============================================================================
# LOAD DATA AND BUILD PREDICTOR LIST
# =============================================================================

message("Loading data: ", data_file)
df         <- readRDS(data_file)
rds_name   <- basename(data_file)

predictors <- feature_set[feature_set %in% names(df)]
if (!length(predictors))
  stop("None of the predictor columns found in the data frame.")

absent <- setdiff(feature_set, names(df))
if (length(absent))
  message("NOTE — absent predictors (will be skipped):\n  ",
          paste(absent, collapse = ", "))

if (MGMT_EVENT != "PI" && !MGMT_EVENT %in% names(df))
  stop("Management column '", MGMT_EVENT, "' not found in ", data_file)

if (PI_ENABLED) {
  pi_cols <- names(df)[grepl("^PI_", names(df))]
  if (!length(pi_cols))
    stop("PI_ENABLED = TRUE but no PI_{gap_label} columns found in ", data_file,
         "\n  Run Script 09 (phytomass_index.R) first.")
  message("PI columns available: ", length(pi_cols))
}


# =============================================================================
# OUTPUT DIRECTORY
# =============================================================================

RESULTS_DIR <- here::here("results", SITE_NAME, "management_effect",
                          MODEL_CHOICE, MGMT_EVENT)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# LOG RUN CONFIGURATION
# =============================================================================

message("\n--- Run configuration ---")
message("  Site          : ", SITE_NAME)
message("  Model         : ", MODEL_CHOICE)
message("  Event         : ", MGMT_EVENT)
message("  PI_ENABLED    : ", PI_ENABLED)
message("  Predictors    : ", paste(predictors, collapse = ", "))
message("  Data rows     : ", nrow(df))
message("  Results dir   : ", RESULTS_DIR)
message("  Model script  : ", model_script)
message("-------------------------\n")


# =============================================================================
# SOURCE THE MODEL SCRIPT
# =============================================================================

run_env             <- new.env(parent = globalenv())
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

# ====================== end run_management_effect.R ==========================