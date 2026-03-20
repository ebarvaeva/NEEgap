# =============================================================================
# run_MDS_miniRECgap.R — Gap-Filling Entry Point: MDS | miniRECgap
# =============================================================================
#
# PURPOSE
#   Entry point for running either of two process-informed gap-filling
#   approaches on a half-hourly eddy-covariance flux dataset:
#
#   MDS — Marginal Distribution Sampling
#     A look-up table method implemented in the REddyProc package
#     (Wutzler et al., 2018).  Fills gaps by finding half-hours with
#     similar meteorological conditions in a moving window around the
#     missing observation.  Suitable for unmanaged sites; no management
#     variables are used.
#
#   miniRECgap — Process-based Reco/GPP Partitioning
#     Fits Lloyd-Taylor Reco and a non-rectangular light response for
#     GPP per regrowth period, then uses these fitted curves to fill gaps.
#     If Grazing_days_since is present in the data, regrowth periods are
#     delineated by grazing events; otherwise a single global period is used.
#     Can be applied to both managed and unmanaged sites.
#
#   Both models produce a single standardised output file
#   (df_cv_all_predictions.rds) compatible with the metrics and plotting
#   scripts used for the machine-learning models.
#
# INPUTS
#   data/{SITE_NAME}.rds
#     Required columns for both models:
#       timestamp  — POSIXct, UTC, 30-min intervals
#       NEE_orig   — measured NEE (µmol m⁻² s⁻¹); NA where gapped
#       PPFD       — photon flux density
#       Temp       — air temperature (°C)
#     Additional columns used when present:
#       VPD               — vapour pressure deficit (MDS)
#       Grazing_days_since — enables per-regrowth partitioning (miniRECgap)
#
# OUTPUTS  (written to results/{SITE_NAME}/{MODEL_CHOICE}/)
#   df_cv_all_predictions.rds
#     The original data frame augmented with prediction columns for each
#     target variable (NEE, Reco, GPP) and gap size (S, M, L, VL).
#
# SECTIONS
#   1. User Settings     — site name, model choice, file paths
#   2. Load Libraries    — REddyProc loaded only when MDS is selected
#   3. Validate          — checks that script and data files exist
#   4. Load Data         — reads the prepared RDS
#   5. Create Output Dir — creates the results folder and logs configuration
#   6. Source Model      — executes the selected model CV script
#
# HOW TO RUN
#   1. Open this script with the project root as the working directory.
#   2. Edit the USER SETTINGS block (Section 1).
#   3. Source the script:  source("scripts/run_MDS_miniRECgap.R")
#
# ADAPTING FOR JC2
#   Change SITE_NAME to "JC2" and update DATA_FILE accordingly.
#   All output directories are derived automatically.
#
# =============================================================================


# =============================================================================
# SECTION 1 — User Settings
# =============================================================================
# Edit the variables below before running.  All other sections are automatic.

# --- Site --------------------------------------------------------------------
SITE_NAME <- "JC1"      # Site identifier; change to "JC2" for the second site

# --- Model choice ------------------------------------------------------------
MODEL_CHOICE <- "miniRECgap"   # "miniRECgap"  or  "MDS"

# --- Input data file ---------------------------------------------------------
# The prepared .rds must contain: timestamp, NEE_orig, PPFD, Temp.
# For miniRECgap, Grazing_days_since is strongly recommended to enable
# per-regrowth-period parameter estimation.
# For MDS, VPD improves look-up accuracy when available.
DATA_FILE <- here::here("data", paste0(SITE_NAME, ".rds"))

# --- Output directory --------------------------------------------------------
# Results are written to:  results/{SITE_NAME}/{MODEL_CHOICE}/
RESULTS_BASE_DIR <- here::here("results")
RESULTS_DIR      <- file.path(RESULTS_BASE_DIR, SITE_NAME, MODEL_CHOICE)

# --- Model script path -------------------------------------------------------
MODEL_SCRIPTS_DIR <- here::here("models")

# =============================================================================
# END OF USER SETTINGS — do not edit below this line for routine use
# =============================================================================


# =============================================================================
# SECTION 2 — Load Libraries
# =============================================================================
# REddyProc is loaded only when MDS is selected, as it has additional system
# dependencies that may not be available on all machines.

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


# =============================================================================
# SECTION 3 — Validate
# =============================================================================
# Confirms that the model script and input data file exist before proceeding.
# Errors at this stage are informative and require no debugging of the model.

if (!MODEL_CHOICE %in% c("miniRECgap", "MDS"))
  stop("MODEL_CHOICE must be 'miniRECgap' or 'MDS'.")

model_script <- file.path(MODEL_SCRIPTS_DIR, paste0(MODEL_CHOICE, "_CV.R"))
if (!file.exists(model_script))
  stop("Model script not found: ", model_script)
if (!file.exists(DATA_FILE))
  stop("Data file not found: ", DATA_FILE)


# =============================================================================
# SECTION 4 — Load Data
# =============================================================================
# The prepared RDS is read and its dimensions printed for confirmation.

message("Loading data: ", DATA_FILE)
df       <- readRDS(DATA_FILE)
rds_name <- basename(DATA_FILE)

message("Data loaded: ", nrow(df), " rows, ", ncol(df), " columns.")
message("Columns: ", paste(names(df), collapse = ", "))


# =============================================================================
# SECTION 5 — Create Output Directory
# =============================================================================
# The results directory is created if it does not already exist.
# Run metadata is printed to the console for verification.

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

message("\n--- Run configuration ---")
message("  Site             : ", SITE_NAME)
message("  Model            : ", MODEL_CHOICE)
message("  Data file        : ", DATA_FILE)
message("  Results dir      : ", RESULTS_DIR)
message("  Model script     : ", model_script)
message("-------------------------\n")


# =============================================================================
# SECTION 6 — Source the Model Script
# =============================================================================
# The model CV script is sourced in an isolated environment to prevent
# variable contamination in the global workspace.  Three variables are
# passed to the model script:
#
#   df          — the full prepared data frame
#   RESULTS_DIR — where all outputs will be written
#   rds_name    — the source filename (logged in run_info.txt)

run_env <- new.env(parent = globalenv())
run_env$df          <- df
run_env$RESULTS_DIR <- RESULTS_DIR
run_env$rds_name    <- rds_name

run_status <- tryCatch({
  source(model_script, local = run_env)
  "completed successfully"
}, error = function(e) {
  paste("ERROR —", conditionMessage(e))
})

message("\n=== Run status: ", run_status, " ===")
message("Results in: ", RESULTS_DIR)

# ========================= end run_MDS_miniRECgap.R ===========================
