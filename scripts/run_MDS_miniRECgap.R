# =============================================================================
# RUN SCRIPT — MDS and miniRECgap Gap-Filling
# =============================================================================
#
# PURPOSE
#   Entry point for running either:
#     MDS        — Marginal Distribution Sampling (REddyProc)  [unmanaged only]
#     miniRECgap — Process-based Reco/GPP partitioning model   [managed or unmanaged]
#
#   Both models produce a single output file:
#     results/{SITE_NAME}/{MODEL_CHOICE}/df_cv_all_predictions.rds
#
# NOTE: MDS and miniRECgap are run on UNMANAGED sites by convention in this
#   project.  miniRECgap can also run on managed sites if Grazing_days_since
#   is present (it then fits per regrowth period rather than globally).
#
# HOW TO ADAPT FOR JC2
#   Change SITE_NAME to "JC2" and update DATA_FILE accordingly.
#   All output directories are derived automatically.
#
# =============================================================================


# =============================================================================
# USER SETTINGS — edit this block for each site / run
# =============================================================================

# --- Site -------------------------------------------------------------------
SITE_NAME <- "JC1"      # change to "JC2" for the second site

# --- Model choice -----------------------------------------------------------
MODEL_CHOICE <- "miniRECgap"   # "miniRECgap"  or  "MDS"

# --- Input data file --------------------------------------------------------
# The .rds must contain: timestamp, NEE_orig, PPFD, Temp
# For miniRECgap: Grazing_days_since is strongly recommended (regrowth periods)
# For MDS:        VPD is used if present; PPFD is used as shortwave radiation (Rg)
DATA_FILE <- here::here("data", paste0(SITE_NAME, ".rds"))
# Adjust the filename above to match your actual file naming convention.

# --- Output directory -------------------------------------------------------
#   Results will be written to:  results/{SITE_NAME}/{MODEL_CHOICE}/
RESULTS_BASE_DIR <- here::here("results")
RESULTS_DIR      <- file.path(RESULTS_BASE_DIR, SITE_NAME, MODEL_CHOICE)

# --- Model script path ------------------------------------------------------
MODEL_SCRIPTS_DIR <- here::here("models")


# =============================================================================
# END OF USER SETTINGS
# =============================================================================


# =============================================================================
# LOAD LIBRARIES
# =============================================================================

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
# VALIDATE
# =============================================================================

if (!MODEL_CHOICE %in% c("miniRECgap", "MDS"))
  stop("MODEL_CHOICE must be 'miniRECgap' or 'MDS'.")

model_script <- file.path(MODEL_SCRIPTS_DIR, paste0(MODEL_CHOICE, "_CV.R"))
if (!file.exists(model_script))
  stop("Model script not found: ", model_script)
if (!file.exists(DATA_FILE))
  stop("Data file not found: ", DATA_FILE)


# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading data: ", DATA_FILE)
df       <- readRDS(DATA_FILE)
rds_name <- basename(DATA_FILE)

message("Data loaded: ", nrow(df), " rows, ", ncol(df), " columns.")
message("Columns: ", paste(names(df), collapse = ", "))


# =============================================================================
# CREATE OUTPUT DIRECTORY
# =============================================================================

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

message("\n--- Run configuration ---")
message("  Site             : ", SITE_NAME)
message("  Model            : ", MODEL_CHOICE)
message("  Data file        : ", DATA_FILE)
message("  Results dir      : ", RESULTS_DIR)
message("  Model script     : ", model_script)
message("-------------------------\n")


# =============================================================================
# SOURCE MODEL SCRIPT
# =============================================================================
# Variables df, RESULTS_DIR, rds_name are passed to the model script.

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
