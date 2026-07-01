# generate_management_csvs.R — extract hard-coded management events into CSVs
#
# Reads the management event tables defined in ManagementEvents.R and writes the
# CSV files the data preparation pipeline expects: grazing and fertiliser event
# windows per site (JC1 2020/2023/2024, JC2 2023/2024), the last event before
# each modelled year, and per-event nitrogen amounts (slurry rows flagged and set
# to 0 N). Run once before run_data_preparation.R; the CSVs can be edited by hand
# afterwards if events change.
#
# Input  : data_preparation/ManagementEvents.R
# Output : data/management_data/*.csv

library(here)
library(dplyr)
library(readr)
library(stringr)
library(lubridate)

# Create output directory for the management CSVs
# Makes data/management_data/ (and parents) only if it does not already exist.
OUT_DIR <- here::here("data", "management_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load the hard-coded management event tables
# Sources ManagementEvents.R, re-trying an explicit path if jc1_events_2020 is missing.
source(here::here("data_preparation", "ManagementEvents.R"), local = TRUE)
if (!exists("jc1_events_2020")) {
  mgmt_path <- here::here("data_preparation", "ManagementEvents.R")
  if (!file.exists(mgmt_path))
    stop("Cannot find ManagementEvents.R. Set the correct path.")
  source(mgmt_path, local = TRUE)
}

# Helper: write one event type (Grazing / Fertiliser) to its own CSV
# Filters the matching events, extracts date / start / end / event_id, sorts, and writes the file.
write_event_csv <- function(events_list, site, event_type, out_dir) {
  combined <- bind_rows(events_list) %>%
    filter(str_detect(event, paste0("^", event_type))) %>%
    mutate(date     = as.character(as.Date(date)),
           event_id = as.integer(str_extract(event, "\\d+$"))) %>%
    select(date, start_time, end_time, event_id) %>%
    arrange(date, start_time)
  
  out_path <- file.path(out_dir,
                        paste0(tolower(event_type), "_events_", site, ".csv"))
  write_csv(combined, out_path)
  invisible(combined)
}

# Write JC1 grazing and fertiliser event CSVs
# Combines the 2020, 2023 and 2024 JC1 event tables into one file per event type.
write_event_csv(list(jc1_events_2020, jc1_events_2023, jc1_events_2024),
                "JC1", "Grazing",    OUT_DIR)
write_event_csv(list(jc1_events_2020, jc1_events_2023, jc1_events_2024),
                "JC1", "Fertiliser", OUT_DIR)

# Write JC2 grazing and fertiliser event CSVs
# Combines the 2023 and 2024 JC2 event tables into one file per event type.
write_event_csv(list(jc2_events_2023, jc2_events_2024), "JC2", "Grazing",    OUT_DIR)
write_event_csv(list(jc2_events_2023, jc2_events_2024), "JC2", "Fertiliser", OUT_DIR)

# JC1 last grazing / fertiliser event before each modelled year
# End time of the final prior-year event, used to seed days-since-event; 2019 values are fallbacks.
last_events_jc1 <- tribble(
  ~year_for, ~event_type,    ~last_end,
  2020L,     "Grazing",      "2019-09-05 00:00:00",   
  2020L,     "Fertiliser",   "2019-09-11 11:00:00",  
  2023L,     "Grazing",      "2022-09-24 00:00:00",
  2023L,     "Fertiliser",   "2022-10-18 00:00:00",
  2024L,     "Grazing",      "2023-10-16 15:00:00",
  2024L,     "Fertiliser",   "2023-05-24 11:00:00"
)
write_csv(last_events_jc1, file.path(OUT_DIR, "last_events_before_year_JC1.csv"))

# JC2 last grazing / fertiliser event before each modelled year
# End time of the final prior-year event, used to seed days-since-event.
last_events_jc2 <- tribble(
  ~year_for, ~event_type,    ~last_end,
  2023L,     "Grazing",      "2022-09-24 00:00:00",
  2023L,     "Fertiliser",   "2022-09-09 00:00:00",
  2024L,     "Grazing",      "2023-10-04 15:00:00",
  2024L,     "Fertiliser",   "2023-09-08 11:00:00"
)
write_csv(last_events_jc2, file.path(OUT_DIR, "last_events_before_year_JC2.csv"))

# JC1 nitrogen applied per fertiliser event (2020 / 2023 / 2024)
# N in kg/ha for each date; slurry rows are flagged is_slurry = TRUE and set to 0 N.
nitrogen_jc1 <- tribble(
  ~fert_date,    ~N_kg_ha, ~is_slurry,
  # ---- 2020 (CAN only; 03/03/2020 slurry excluded) ----
  "02/04/2020",   50.0,    FALSE,
  "11/05/2020",   40.0,    FALSE,
  "03/06/2020",   27.0,    FALSE,
  "29/06/2020",   20.0,    FALSE,
  "14/08/2020",   27.0,    FALSE,
  "14/09/2020",   27.0,    FALSE,
  # ---- 2023 ----
  "20/02/2023",   18.5,    FALSE,
  "19/05/2023",    0.0,    TRUE,    # Slurry — excluded from N predictor
  # ---- 2024 ----
  "02/04/2024",   30.4,    FALSE,
  "07/05/2024",   23.0,    FALSE,
  "29/05/2024",   17.0,    FALSE,
  "20/06/2024",   10.0,    FALSE,
  "06/08/2024",   10.0,    FALSE,
  "22/08/2024",    0.0,    TRUE     # Slurry — excluded from N predictor
)
write_csv(nitrogen_jc1, file.path(OUT_DIR, "nitrogen_amounts_JC1.csv"))

# JC2 nitrogen applied per fertiliser event (2023 / 2024)
# N in kg/ha for each date; slurry rows are flagged is_slurry = TRUE and set to 0 N.
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
  "19/04/2024",    0.0,    TRUE,    # Slurry — excluded from N predictor
  "07/05/2024",   30.4,    FALSE,
  "29/05/2024",   25.0,    FALSE,
  "19/06/2024",   28.0,    FALSE,
  "10/07/2024",   20.0,    FALSE,
  "06/08/2024",   25.0,    FALSE
)
write_csv(nitrogen_jc2, file.path(OUT_DIR, "nitrogen_amounts_JC2.csv"))
