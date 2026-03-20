# =============================================================================
# generate_management_csvs.R — ONE-TIME SETUP SCRIPT
# =============================================================================
# PURPOSE
#   Extract management event data (hard-coded in ManagementEvents.R) into
#   the CSV files that the data preparation pipeline expects.
#
#   Run this script ONCE before running run_data_preparation.R.
#   After running it, you can edit the CSVs manually if events change.
#
# OUTPUTS  (all written to data/management_data/)
#   grazing_events_JC1.csv        — JC1 grazing windows 2023–2024
#   grazing_events_JC2.csv        — JC2 grazing windows 2023–2024
#   fertiliser_events_JC1.csv     — JC1 fertiliser windows 2023–2024
#   fertiliser_events_JC2.csv     — JC2 fertiliser windows 2023–2024
#   last_events_before_year_JC1.csv — last grazing/fertiliser end before 2023 and 2024
#   last_events_before_year_JC2.csv
#   nitrogen_amounts_JC1.csv      — N kg ha⁻¹ per mineral fertiliser event
#   nitrogen_amounts_JC2.csv
# =============================================================================

library(here)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)

OUT_DIR <- here::here("data", "management_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source ManagementEvents.R to load jc1/jc2 event tables
# Adjust path if the file is stored elsewhere
source(here::here("data_preparation", "ManagementEvents.R"),
       local = TRUE)
# Fallback: if above fails, try direct path
if (!exists("jc1_events_2023")) {
  mgmt_path <- here::here("data_preparation", "ManagementEvents.R")
  if (!file.exists(mgmt_path))
    stop("Cannot find ManagementEvents.R. Set the correct path below:\n",
         "  source(<path to ManagementEvents.R>, local = TRUE)")
  source(mgmt_path, local = TRUE)
}

# =============================================================================
# HELPER: split event table into separate event-type CSVs
# =============================================================================
write_event_csv <- function(events_list, site, event_type, out_dir) {
  # Bind all years, filter to the event type, keep columns: date, start_time, end_time
  combined <- bind_rows(events_list) %>%
    filter(str_detect(event, paste0("^", event_type))) %>%
    mutate(date = as.character(as.Date(date)),
           event_id = as.integer(str_extract(event, "\\d+$"))) %>%
    select(date, start_time, end_time, event_id) %>%
    arrange(date, start_time)

  out_path <- file.path(out_dir, paste0(tolower(event_type), "_events_", site, ".csv"))
  write_csv(combined, out_path)
  message("  Written: ", basename(out_path), "  (", nrow(combined), " rows)")
  invisible(combined)
}

# =============================================================================
# JC1
# =============================================================================
message("Writing JC1 management event CSVs ...")
write_event_csv(list(jc1_events_2023, jc1_events_2024), "JC1", "Grazing",    OUT_DIR)
write_event_csv(list(jc1_events_2023, jc1_events_2024), "JC1", "Fertiliser", OUT_DIR)

# =============================================================================
# JC2
# =============================================================================
message("Writing JC2 management event CSVs ...")
write_event_csv(list(jc2_events_2023, jc2_events_2024), "JC2", "Grazing",    OUT_DIR)
write_event_csv(list(jc2_events_2023, jc2_events_2024), "JC2", "Fertiliser", OUT_DIR)

# =============================================================================
# Last events before each year  (hardcoded from CreateData / Corrected_timestamps)
# =============================================================================
# These values come from:
#   01. Corrected_timestamps_JC1.R  lines 579-583
#   01. CreateData_JC1_2023.R       lines 819-831
#   01. Corrected_timestamps_JC2.R  lines 606-610
# Edit here if your last pre-year events differ.

message("Writing last_events_before_year CSVs ...")

last_events_jc1 <- tribble(
  ~year_for, ~event_type,    ~last_end,
  2023L,     "Grazing",      "2022-09-24 00:00:00",
  2023L,     "Fertiliser",   "2022-10-18 00:00:00",
  2024L,     "Grazing",      "2023-10-16 15:00:00",
  2024L,     "Fertiliser",   "2023-05-24 11:00:00"
)
write_csv(last_events_jc1, file.path(OUT_DIR, "last_events_before_year_JC1.csv"))
message("  Written: last_events_before_year_JC1.csv")

last_events_jc2 <- tribble(
  ~year_for, ~event_type,    ~last_end,
  2023L,     "Grazing",      "2022-09-24 00:00:00",
  2023L,     "Fertiliser",   "2022-09-09 00:00:00",
  2024L,     "Grazing",      "2023-10-04 15:00:00",
  2024L,     "Fertiliser",   "2023-09-08 11:00:00"
)
write_csv(last_events_jc2, file.path(OUT_DIR, "last_events_before_year_JC2.csv"))
message("  Written: last_events_before_year_JC2.csv")

# =============================================================================
# Nitrogen amounts  (hardcoded from 01e. Regrowth_period.R  N_windows_jc1/jc2)
# Slurry events are marked is_slurry = TRUE and excluded from the N predictor.
# Edit this section if your nitrogen schedule changes.
# =============================================================================
message("Writing nitrogen_amounts CSVs ...")

nitrogen_jc1 <- tribble(
  ~fert_date,    ~N_kg_ha, ~is_slurry,
  "20/02/2023",   18.5,    FALSE,
  "19/05/2023",    0.0,    TRUE,   # Slurry — excluded from N predictor
  "02/04/2024",   30.4,    FALSE,
  "07/05/2024",   23.0,    FALSE,
  "29/05/2024",   17.0,    FALSE,
  "20/06/2024",   10.0,    FALSE,
  "06/08/2024",   10.0,    FALSE,
  "22/08/2024",    0.0,    TRUE    # Slurry — excluded from N predictor
)
write_csv(nitrogen_jc1, file.path(OUT_DIR, "nitrogen_amounts_JC1.csv"))
message("  Written: nitrogen_amounts_JC1.csv")

nitrogen_jc2 <- tribble(
  ~fert_date,    ~N_kg_ha, ~is_slurry,
  "20/02/2023",   30.0,    FALSE,
  "04/04/2023",   35.0,    FALSE,
  "26/04/2023",   32.0,    FALSE,
  "22/05/2023",   30.0,    FALSE,
  "20/06/2023",   30.0,    FALSE,
  "12/07/2023",   27.0,    FALSE,
  "24/08/2023",   27.0,    FALSE,
  "08/09/2023",   21.0,    FALSE,
  "02/04/2024",   38.0,    FALSE,
  "19/04/2024",    0.0,    TRUE,   # Slurry — excluded from N predictor
  "07/05/2024",   30.4,    FALSE,
  "29/05/2024",   25.0,    FALSE,
  "19/06/2024",   28.0,    FALSE,
  "10/07/2024",   20.0,    FALSE,
  "06/08/2024",   25.0,    FALSE
)
write_csv(nitrogen_jc2, file.path(OUT_DIR, "nitrogen_amounts_JC2.csv"))
message("  Written: nitrogen_amounts_JC2.csv")

# =============================================================================
# SUMMARY
# =============================================================================
message("\nAll management CSVs written to: ", OUT_DIR)
message("Files created:")
for (f in list.files(OUT_DIR, pattern = "\\.csv$")) message("  ", f)
message("\nVerify dates and amounts before running run_data_preparation.R.")
