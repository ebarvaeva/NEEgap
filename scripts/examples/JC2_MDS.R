# JC2_MDS.R — Run script for MDS NEE gap-filling
#
# Entry point for one site run. It loads the prepared dataset, creates the
# output directory, and sources models/{MODEL_CHOICE}_CV.R to gap-fill across
# four gap-size classes (VL = very large, L = large, M = medium, S = small):
#   MDS        — Marginal Distribution Sampling (REddyProc); fills NEE only
#   miniRECgap — process model NEE = Reco - GPP; fills NEE, Reco and GPP
#
# Required input columns (in DATA_FILE):
#   timestamp   — POSIXct / parseable half-hourly date-time
#   NEE_orig    — net ecosystem exchange (µmol m-2 s-1)
#   PPFD, Temp  — light and air-temperature drivers
#   miniRECgap also uses Grazing_days_since (managed sites) to fit per regrowth
#   period; MDS uses VPD if present and treats PPFD as shortwave radiation (Rg).
#
# Output (in RESULTS_DIR):
#   df_cv_all_predictions.rds   — gap-filled predictions (4 columns: VL/L/M/S)
#
# To adapt: change SITE_NAME / MODEL_CHOICE and point DATA_FILE at the new .rds.


# Edit these settings for each site / run:
SITE_NAME    <- "JC2"        # site identifier ("JC2", "JC3", …)
MODEL_CHOICE <- "MDS"   # "MDS" or "miniRECgap"; names the model script and results subfolder

# prepared cross-validation dataset
DATA_FILE <- here::here("data/data_prepared", paste0(SITE_NAME, "_cv.rds"))

# output directory: results/{SITE}/{MODEL}
RESULTS_DIR <- here::here("results", SITE_NAME, MODEL_CHOICE)

# directory holding the model scripts
MODEL_SCRIPTS_DIR <- here::here("models")


# attach the packages used across both scripts (REddyProc is only needed for MDS)
suppressPackageStartupMessages({
  library(here); library(dplyr); library(lubridate); library(tibble)
})
if (MODEL_CHOICE == "MDS") library(REddyProc)

# check the model choice is valid and both inputs exist
if (!MODEL_CHOICE %in% c("MDS", "miniRECgap")) stop("MODEL_CHOICE must be 'MDS' or 'miniRECgap'.")
model_script <- file.path(MODEL_SCRIPTS_DIR, paste0(MODEL_CHOICE, "_CV.R"))
if (!file.exists(model_script)) stop("Model script not found: ", model_script)
if (!file.exists(DATA_FILE))    stop("Data file not found: ", DATA_FILE)

# load the prepared data
df <- readRDS(DATA_FILE)

# create the results directory (no-op if it already exists)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# pass settings to the model script via its own environment, then source it
run_env <- new.env(parent = globalenv())
run_env$df          <- df
run_env$RESULTS_DIR <- RESULTS_DIR

source(model_script, local = run_env)
