# 01_run_data_preparation.R — master runner for the data-preparation pipeline
#
# Sources sub-scripts 02-07 in order for each site x year, injecting SITE, YEAR
# and flags into each script's own environment. JC1 has an extra 2020 open-path
# year that is handled before the standard 2023/2024 loop. Scripts 04-06
# (management) run only when HAS_MANAGEMENT is TRUE for the site. This runner
# writes no output itself; each sub-script saves its own results.

library(here)
source("data_preparation/generate_management_csvs.R")

# sites and standard years to process
SITES <- c("JC1", "JC2")
YEARS <- c(2023, 2024)

# per-site management switch (FALSE skips scripts 04-06)
HAS_MANAGEMENT <- list(JC1 = TRUE, JC2 = TRUE)

# sub-scripts sourced by this runner
SCRIPT_DIR <- here::here("data_preparation")
scripts <- list(
  `02` = file.path(SCRIPT_DIR, "02_timestamp_correction.R"),
  `03` = file.path(SCRIPT_DIR, "03_predictors_meteo.R"),
  `04` = file.path(SCRIPT_DIR, "04_management_days_since.R"),
  `05` = file.path(SCRIPT_DIR, "05_nitrogen.R"),
  `06` = file.path(SCRIPT_DIR, "06_grass_canopy.R"),
  `07` = file.path(SCRIPT_DIR, "07_bind_years.R")
)

# stop early if any sub-script file is missing
missing_scripts <- names(scripts)[!vapply(scripts, file.exists, logical(1))]
if (length(missing_scripts))
  stop("Scripts not found: ", paste(missing_scripts, collapse=", "))

# source one sub-script with SITE, YEAR and flags injected into its environment
run_script <- function(script_path, site, year, twin_site = NULL,
                       has_mgmt = TRUE, is_open_path = FALSE) {
  env <- new.env(parent = globalenv())
  env$SITE          <- site
  env$YEAR          <- year
  env$HAS_MANAGEMENT <- has_mgmt
  env$IS_OPEN_PATH  <- is_open_path   # tells script 03 to use the 2020 branch
  if (!is.null(twin_site)) env$TWIN_SITE <- twin_site
  source(script_path, local = env)
  invisible(NULL)
}

# JC1 2020 open-path year — run before the standard loop, only if its QC exists
{
  site <- "JC1"
  yr   <- 2020
  twin <- "JC2"
  mgmt <- isTRUE(HAS_MANAGEMENT[[site]])
  
  qc_input <- here::here("data", "data_qc",
                         paste0("results_qc_", site, "_", yr), "merged_qc.rds")
  if (!file.exists(qc_input)) {
    # QC output for JC1 2020 not found — skip this year
  } else {
    tryCatch({
      run_script(scripts[["02"]], site = site, year = yr)
      run_script(scripts[["03"]], site = site, year = yr,
                 twin_site = twin, is_open_path = TRUE)   # open-path flag
      if (mgmt) {
        run_script(scripts[["04"]], site = site, year = yr)
        run_script(scripts[["05"]], site = site, year = yr)
        run_script(scripts[["06"]], site = site, year = yr)
      }
    }, error = function(e) {
      # a failure here skips JC1 2020; the rest of the pipeline still runs
    })
  }
}

# standard years (2023 + 2024) for every site
for (site in SITES) {
  twin <- setdiff(SITES, site)
  mgmt <- isTRUE(HAS_MANAGEMENT[[site]])
  
  for (yr in YEARS) {
    qc_input <- here::here("data", "data_qc",
                           paste0("results_qc_", site, "_", yr), "merged_qc.csv")
    # skip this site/year if its QC output is missing
    if (!file.exists(qc_input)) next
    
    tryCatch({
      run_script(scripts[["02"]], site = site, year = yr)
      run_script(scripts[["03"]], site = site, year = yr, twin_site = twin)
      
      if (mgmt) {
        run_script(scripts[["04"]], site = site, year = yr)
        run_script(scripts[["05"]], site = site, year = yr)
        run_script(scripts[["06"]], site = site, year = yr)
      }
    }, error = function(e) {
      # a failing year is skipped so the other years still run
    })
  }
  
  # bind all available years for this site (includes 2020 for JC1)
  run_script(scripts[["07"]], site = site, year = NA, has_mgmt = mgmt)
}
