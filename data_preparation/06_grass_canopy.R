# =============================================================================
# 06_grass_canopy.R -- Grass Height and Dry Biomass
# =============================================================================
#
# PURPOSE
#   Creates two canopy structure predictors at half-hourly resolution:
#     grass_height   (cm)          -- sward height
#     grass_biomass  (kg DM ha-1)  -- above-ground dry matter
#
# WHY CANOPY STRUCTURE MATTERS
#   Leaf area index (LAI) and standing biomass are the primary determinants
#   of canopy light interception and hence gross primary production (GPP).
#   In managed grassland, canopy structure is disrupted by grazing and
#   recovers over days to weeks.  Providing the model with grass_height and
#   grass_biomass as predictors allows it to condition GPP and respiration
#   estimates on the current canopy state, which is not captured by radiation
#   or temperature alone.
#
# INTERPOLATION DESIGN
#   Field measurements are taken weekly or bi-weekly; between measurement
#   dates the true canopy state is unknown.  Linear interpolation (zoo::
#   na.approx, rule = 2) is used to fill the gaps.  Before interpolation,
#   all grazing days are anchored to POST_GRAZING_HEIGHT_CM (default 4 cm),
#   the observed post-grazing residual.  This anchor prevents the interpolation
#   from extrapolating a smooth curve through a grazing event, which would
#   underestimate the abrupt structural change caused by defoliation.
#   rule = 2 (constant extrapolation at boundaries) fills the start and end
#   of the year with the nearest observed value rather than producing NA.
#
# INPUTS
#   data/data_prepared/{SITE}_{YEAR}_predictors.rds
#     Must already have Grazing_days_since (Script 04) so that grazing dates
#     can be identified as dates where Grazing_days_since == 0.
#
#   data/canopy_data/canopy_measurements_{SITE}_{YEAR}.csv
#     One row per measurement date × subplot.  Columns:
#       date          — YYYY-MM-DD
#       subplot_id    — subplot/pasture identifier (character or integer)
#       height_cm     — grass height (cm)
#       drymass_kg_ha — dry biomass (kg DM ha⁻¹)
#
# OUTPUTS
#   data/data_prepared/{SITE}_{YEAR}_predictors.rds  (updated in-place with
#     grass_height and grass_biomass)
# =============================================================================

library(here)
library(dplyr)
library(readr)
library(lubridate)
library(tidyr)
library(zoo)

# -----------------------------------------------------------------------------
# USER SETTINGS
# -----------------------------------------------------------------------------
TZ   <- "UTC"

# Height assigned to all grazing days before interpolation.
# This value corresponds to the typical post-grazing residual height observed
# at the JC sites.  Adjust if measurements at a new site differ.
POST_GRAZING_HEIGHT_CM <- 4   # cm

# -----------------------------------------------------------------------------
# SECTION 1 -- Load predictors RDS and identify grazing dates
# -----------------------------------------------------------------------------
# Loads the per-year predictors RDS (updated by Scripts 03-05) and identifies
# all dates on which grazing occurred.  Grazing dates are read directly from
# the event CSV rather than inferred from Grazing_days_since == 0, because
# the CSV is the authoritative source and covers all years (including the
# previous year, which contributes the first-of-year days_since value).
# The resulting grazing_dates vector is used in Section 3 to set the canopy
# anchor before interpolation.
in_path <- here::here("data", "data_prepared",
                       paste0(SITE, "_", YEAR, "_predictors.rds"))
if (!file.exists(in_path)) stop("Predictors RDS not found: ", in_path)
df <- readRDS(in_path) %>%
  mutate(timestamp = as.POSIXct(timestamp, tz = TZ)) %>%
  arrange(timestamp)
message("Loaded: ", in_path)

# Derive grazing dates from the events CSV (ALL years, not just current year)
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
message("  Grazing dates identified (all years): ", length(grazing_dates))

# -----------------------------------------------------------------------------
# SECTION 2 -- Load canopy measurements
# -----------------------------------------------------------------------------
# Reads the shared canopy RDS (all sites x years in wide format) and selects
# the height and drymass columns for the current site.  Measured heights below
# POST_GRAZING_HEIGHT_CM are clamped to that minimum to prevent the interpolation
# from dipping below the post-grazing floor between measurement dates.
canopy_path <- here::here("data", "canopy_data",
                            "grass_drybiomass_JC_2023_2024.rds")
if (!file.exists(canopy_path)) stop("Canopy RDS not found: ", canopy_path)

canopy_all <- readRDS(canopy_path) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

# Select columns for this site
height_col  <- paste0("height_",  SITE)
drymass_col <- paste0("drymass_", SITE)
for (col in c(height_col, drymass_col)) {
  if (!col %in% names(canopy_all))
    stop("Column '", col, "' not found in canopy RDS. ",
         "Available: ", paste(names(canopy_all), collapse = ", "))
}

canopy_daily <- canopy_all %>%
  transmute(
    date          = date,
    height_cm     = .data[[height_col]],
    drymass_kg_ha = .data[[drymass_col]]
  ) %>%
  mutate(
    # Clamp any measured height below minimum (matches old pipeline)
    height_cm = if_else(!is.na(height_cm) & height_cm < POST_GRAZING_HEIGHT_CM,
                        as.numeric(POST_GRAZING_HEIGHT_CM), height_cm)
  )

message("  Canopy measurement dates for ", SITE, " ", YEAR, ": ",
        sum(!is.na(canopy_daily$height_cm)))

# -----------------------------------------------------------------------------
# SECTION 3 -- Build full daily grid and apply grazing anchor
# -----------------------------------------------------------------------------
# A full daily grid spanning all years covered by the canopy data is constructed
# and left-joined with the measurement data.  The grazing anchor (POST_GRAZING_
# HEIGHT_CM) is then applied to every date in grazing_dates, overwriting both
# measured and missing-but-would-be-interpolated values on those days.
# This must happen BEFORE interpolation so that the anchor acts as a knot in
# the piecewise linear interpolation, not just as an initial condition.
daily_grid <- tibble(
  date = seq(
    as.Date(paste0(min(lubridate::year(canopy_daily$date), na.rm = TRUE), "-01-01")),
    as.Date(paste0(max(lubridate::year(canopy_daily$date), na.rm = TRUE), "-12-31")),
    by = "1 day"
  )
)

canopy_full <- daily_grid %>%
  left_join(canopy_daily, by = "date") %>%
  mutate(
    height_cm = if_else(date %in% grazing_dates,
                        as.numeric(POST_GRAZING_HEIGHT_CM), height_cm)
  ) %>%
  arrange(date)

# -----------------------------------------------------------------------------
# SECTION 4 -- Linear interpolation across the daily grid
# -----------------------------------------------------------------------------
# zoo::na.approx() interpolates linearly between non-NA anchor points.
# Interpolation is performed within each calendar year separately (group_by year)
# so that the year-end value of 2023 does not influence the year-start of 2024.
# rule = 2 extends the boundary values (constant extrapolation) rather than
# leaving NAs at the edges of the year when no measurement exists at Jan 1 or
# Dec 31.  The same logic applies to drymass.
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


message("  NA height after interpolation:  ",
        sum(is.na(canopy_interp$height_cm)))
message("  NA drymass after interpolation: ",
        sum(is.na(canopy_interp$drymass_kg_ha)))

# -----------------------------------------------------------------------------
# SECTION 5 -- Expand daily values to half-hourly resolution and join
# -----------------------------------------------------------------------------
# Canopy state changes on a daily timescale, so the same interpolated value is
# assigned to all 48 half-hours within each calendar day.  The join key is
# as.Date(timestamp), which aligns UTC timestamps to the day they fall in.
df <- df %>%
  mutate(date = as.Date(timestamp)) %>%
  left_join(
    canopy_interp %>% select(date, grass_height = height_cm,
                             grass_biomass = drymass_kg_ha),
    by = "date"
  ) %>%
  select(-date)

message("  grass_height NA after join:  ", sum(is.na(df$grass_height)))
message("  grass_biomass NA after join: ", sum(is.na(df$grass_biomass)))

# -----------------------------------------------------------------------------
# SECTION 6 -- Save updated RDS
# -----------------------------------------------------------------------------
# Overwrites the per-year predictors RDS, adding grass_height and grass_biomass.
# After this script, the RDS contains all BASE and management predictors.
# Script 07 will bind 2023 and 2024 into the final site-level dataset.
saveRDS(df, in_path)
message("Updated RDS saved: ", in_path)
message("\nScript 06 complete.")
