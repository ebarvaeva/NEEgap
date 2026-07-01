# 06_grass_canopy.R — grass_height (cm) & grass_biomass (kg DM ha-1)
#
# Builds two half-hourly canopy predictors from sparse field measurements.
# Grazing days are anchored to a fixed post-grazing height (4 cm) before
# interpolation, sparse daily values are linearly interpolated (per year), and
# the daily series is then expanded to 30-min resolution by a date join.
#
# Inputs :
#   data/data_prepared/{SITE}_{YEAR}_predictors.rds        (must have Grazing_days_since)
#   data/management_data/grazing_events_{SITE}.csv         (grazing dates -> post-grazing anchors)
#   data/canopy_data/grass_drybiomass_JC_2023_2024.rds     (date, height_{SITE}, drymass_{SITE})
# Output :
#   data/data_prepared/{SITE}_{YEAR}_predictors.rds  (updated in-place)

library(here)
library(dplyr)
library(readr)
library(lubridate)
library(tidyr)
library(zoo)

TZ   <- "UTC"

# height (cm) assigned on grazing days before interpolation
POST_GRAZING_HEIGHT_CM <- 4

# load the predictors RDS for this site-year
in_path <- here::here("data", "data_prepared",
                      paste0(SITE, "_", YEAR, "_predictors.rds"))
if (!file.exists(in_path)) stop("Predictors RDS not found: ", in_path)
df <- readRDS(in_path) %>%
  mutate(timestamp = as.POSIXct(timestamp, tz = TZ)) %>%
  arrange(timestamp)

# grazing dates from the events CSV (all years) — used as post-grazing anchors
path_graz_csv <- here::here("data", "management_data",
                            paste0("grazing_events_", SITE, ".csv"))
if (file.exists(path_graz_csv)) {
  graz_csv <- readr::read_csv(path_graz_csv, show_col_types = FALSE) %>%
    mutate(
      date       = as.Date(date),
      start_time = stringr::str_pad(start_time, 5, pad = "0"),
      end_time   = stringr::str_pad(end_time,   5, pad = "0"),
      start_ct   = as.POSIXct(paste(date, start_time), tz = "UTC"),
      end_ct     = dplyr::if_else(
        end_time == "24:00",
        as.POSIXct(paste(date + lubridate::days(1), "00:00"), tz = "UTC"),
        as.POSIXct(paste(date, end_time), tz = "UTC")
      )
    ) %>%
    filter(!is.na(start_ct), end_ct > start_ct)
  
  grazing_dates <- purrr::map_dfr(seq_len(nrow(graz_csv)), function(i) {
    tibble(date = seq.Date(as.Date(graz_csv$start_ct[i]),
                           as.Date(graz_csv$end_ct[i] - lubridate::seconds(1)),
                           by = "1 day"))
  }) %>% distinct() %>% pull(date)
} else {
  warning("grazing_events CSV not found — no post-grazing anchors applied.")
  grazing_dates <- as.Date(character())
}

# canopy measurements: single RDS holding all sites x years
canopy_path <- here::here("data", "canopy_data",
                          "grass_drybiomass_JC_2023_2024.rds")
if (!file.exists(canopy_path)) stop("Canopy RDS not found: ", canopy_path)

canopy_all <- readRDS(canopy_path) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

# select this site's height/biomass columns
height_col  <- paste0("height_",  SITE)
drymass_col <- paste0("drymass_", SITE)
for (col in c(height_col, drymass_col)) {
  if (!col %in% names(canopy_all))
    stop("Column '", col, "' not found in canopy RDS. ",
         "Available: ", paste(names(canopy_all), collapse = ", "))
}

# clamp any measured height below the post-grazing minimum
canopy_daily <- canopy_all %>%
  transmute(
    date          = date,
    height_cm     = .data[[height_col]],
    drymass_kg_ha = .data[[drymass_col]]
  ) %>%
  mutate(
    height_cm = if_else(!is.na(height_cm) & height_cm < POST_GRAZING_HEIGHT_CM,
                        as.numeric(POST_GRAZING_HEIGHT_CM), height_cm)
  )

# full daily grid over the measurement years
daily_grid <- tibble(
  date = seq(
    as.Date(paste0(min(lubridate::year(canopy_daily$date), na.rm = TRUE), "-01-01")),
    as.Date(paste0(max(lubridate::year(canopy_daily$date), na.rm = TRUE), "-12-31")),
    by = "1 day"
  )
)

# apply post-grazing height on grazing days
canopy_full <- daily_grid %>%
  left_join(canopy_daily, by = "date") %>%
  mutate(
    height_cm = if_else(date %in% grazing_dates,
                        as.numeric(POST_GRAZING_HEIGHT_CM), height_cm)
  ) %>%
  arrange(date)

# linear interpolation across the daily grid, within each year
canopy_interp <- canopy_full %>%
  mutate(year = lubridate::year(date)) %>%
  dplyr::group_by(year) %>%
  dplyr::mutate(
    height_cm     = zoo::na.approx(height_cm,     x = as.numeric(date),
                                   na.rm = FALSE, rule = 2),
    drymass_kg_ha = zoo::na.approx(drymass_kg_ha, x = as.numeric(date),
                                   na.rm = FALSE, rule = 2)
  ) %>%
  dplyr::ungroup() %>%
  select(-year)

# expand daily values to half-hourly rows via a date join
df <- df %>%
  mutate(date = as.Date(timestamp)) %>%
  left_join(
    canopy_interp %>% select(date, grass_height = height_cm,
                             grass_biomass = drymass_kg_ha),
    by = "date"
  ) %>%
  select(-date)

# save updated RDS in place
saveRDS(df, in_path)
