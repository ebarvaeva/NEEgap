# run_MDS_miniRECgap.R — Run script for MDS / miniRECgap NEE gap-filling
# Entry point for one site run of a process-based model. It loads the input
# dataset, creates the output directory, and sources the chosen model script
# (models/MDS_CV.R or models/miniRECgap_CV.R) to gap-fill NEE across the four
# gap-size classes (VL = very large, L = large, M = medium, S = small).
#   MDS        — Marginal Distribution Sampling (REddyProc); unmanaged sites,
#                needs a complete annual timeline.
#   miniRECgap — process-based Reco/GPP model; managed or unmanaged. With
#                Grazing_days_since it fits per regrowth period, else globally.
# Required input columns (in DATA_FILE):
#   timestamp — POSIXct / parseable half-hourly date-time
#   NEE_orig  — net ecosystem exchange (µmol m-2 s-1)
#   PPFD      — photon flux density (used as shortwave radiation by MDS)
#   Temp      — air temperature
#   Grazing_days_since — optional; enables per-regrowth fitting in miniRECgap
#   VPD       — optional; used by MDS if present
# Writes: nothing itself — the sourced model script saves predictions to
#   results/{SITE_NAME}/{MODEL_CHOICE}/df_cv_all_predictions.rds
# To adapt: set SITE_NAME and MODEL_CHOICE, and point DATA_FILE at the new .rds.

# User settings — edit for each run
SITE_NAME    <- "JC1"          # "JC1" or "JC2"
MODEL_CHOICE <- "miniRECgap"   # "miniRECgap" or "MDS"

# Input data file (must contain timestamp, NEE_orig, PPFD, Temp)
DATA_FILE <- here::here("data", paste0(SITE_NAME, ".rds"))

# Output directory: results/{SITE_NAME}/{MODEL_CHOICE}/
RESULTS_BASE_DIR <- here::here("results")
RESULTS_DIR      <- file.path(RESULTS_BASE_DIR, SITE_NAME, MODEL_CHOICE)

# Model script location
MODEL_SCRIPTS_DIR <- here::here("models")

# Packages (REddyProc is loaded only for MDS)
suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(lubridate)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(glue)
})

if (MODEL_CHOICE == "MDS") {
  if (!requireNamespace("REddyProc", quietly = TRUE))
    stop("Package 'REddyProc' is required for MDS.  Run: install.packages('REddyProc')")
  library(REddyProc)
}

# Validate settings and locate the model CV script
if (!MODEL_CHOICE %in% c("miniRECgap", "MDS"))
  stop("MODEL_CHOICE must be 'miniRECgap' or 'MDS'.")

model_script <- file.path(MODEL_SCRIPTS_DIR, paste0(MODEL_CHOICE, "_CV.R"))
if (!file.exists(model_script))
  stop("Model script not found: ", model_script)
if (!file.exists(DATA_FILE))
  stop("Data file not found: ", DATA_FILE)

# Load the input data
df       <- readRDS(DATA_FILE)
rds_name <- basename(DATA_FILE)

# Create the output directory
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Pass data + settings to the model script in an isolated environment
run_env <- new.env(parent = globalenv())
run_env$df          <- df
run_env$RESULTS_DIR <- RESULTS_DIR
run_env$rds_name    <- rds_name

# Run the chosen model script (it writes the predictions to RESULTS_DIR)
source(model_script, local = run_env)
