# JC1_2020_qc.R — quality-control the half-hourly CO2 flux for JC1 2020 (open-path)
#
# The 2020 open-path setup provides only eddypro + biomet exports, plus a separate
# canopy/displacement/roughness file (extra) and an AGC file. Aligns them on a
# half-hourly UTC timestamp and merges. Applies an AGC gain filter and a night-time
# uptake filter, the base range / random-error filters, a seasonal u* threshold
# (REddyProc), the 1-9 stationarity flag levels, then the 2-D footprint filter
# (Kljun 2015) against the site boundary polygon at 70/80/90 % thresholds.
#
# Input  : data/raw_data/JC1_2020_{eddypro,biomet,extra,AGC}.xlsx
#          data_preparation/quality_control_nee/grid_qc/nasco_site_list_wkt.csv
# Output : data/data_qc/results_qc_JC1_2020/merged_qc.csv
#          data/data_qc/results_qc_JC1_2020/merged_qc.rds
#----------------------------------------------------------
# Packages
library(here)
library(readr)
library(readxl)
library(tidyverse)
library(hms)
library(REddyProc)
library(zoo)
library(ranger)
#----------------------------------------------------------
# Paths 

eddypro_path <- here::here("data/raw_data/JC1_2020_eddypro.xlsx")
biomet_path  <- here::here("data/raw_data/JC1_2020_biomet.xlsx")
solar_path <- here::here("data", "meteireann_data", "solar_hourly.xlsx")
extra_path <- here::here("data/raw_data/JC1_2020_extra.xlsx")
agc_path<- here::here("data/raw_data/JC1_2020_AGC.xlsx")

# Year
year_no <- 2020

# make sure the folder exists
results_path <- here::here("data/data_qc/results_qc_JC1_2020")

#----------------------------------------------------------
# define thresholds 
th_co2_low      <- -40
th_co2_high     <-  30
th_fetch70      <- NA_real_   # e.g. 100 to require FETCH_70_fluxnet >= 100
th_fetch80      <- NA_real_
th_fetch90      <- NA_real_
th_rand_err     <- 100        # max random error
th_agc          <- 55
# th_sig_strength <- 70         # min signal strength
#th_flow_low     <- 12         # valid flow range low
#th_flow_high    <- 18         # valid flow range high
# th_ppfd_day     <- 10         # ppfd daytime threshold
# th_ustar        <- 0.1

#----------------------------------------------------------
# define polygon name for footprint filter
site_key <- "Johnstown Castle - Teagasc (1)"

#----------------------------------------------------------
# Read data (the if else will choose either read_csv() or read_excel())

## eddypro
if (grepl("\\.csv$", eddypro_path,  ignore.case = TRUE)) {
  eddypro <- readr::read_csv(eddypro_path)
} else if (grepl("\\.xlsx?$", eddypro_path, ignore.case = TRUE)) {
  eddypro <- readxl::read_excel(eddypro_path)
} else stop("Unsupported file type for eddypro_path")

## biomet
if (grepl("\\.csv$", biomet_path,  ignore.case = TRUE)) {
  biomet <- readr::read_csv(biomet_path)
} else if (grepl("\\.xlsx?$", biomet_path, ignore.case = TRUE)) {
  biomet <- readxl::read_excel(biomet_path)
} else stop("Unsupported file type for biomet_path")

#----------------------------------------------------------
# Timestamp

eddypro <- eddypro %>%
  mutate(
    timestamp = update(ymd_hms(time, tz = "UTC"),
                       year  = year(as_date(date)),
                       month = month(as_date(date)),
                       mday  = mday(as_date(date)))
  )  %>%
  relocate(timestamp, .before = 1)

ts <- eddypro$timestamp
for (i in which(is.na(ts))) {
  j <- max(which(!is.na(ts[seq_len(i-1)])), na.rm = TRUE)
  if (is.finite(j)) ts[i] <- ts[j] + minutes(30)   
}
eddypro$timestamp <- ts

biomet <- biomet %>%
  mutate(
    timestamp = update(ymd_hms(time, tz = "UTC"),
                       year  = year(as_date(date)),
                       month = month(as_date(date)),
                       mday  = mday(as_date(date)))
  )  %>%
  relocate(timestamp, .before = 1)

ts <- biomet$timestamp
for (i in seq_along(ts)) if (i > 1 && is.na(ts[i]) && !is.na(ts[i-1])) ts[i] <- ts[i-1] + minutes(30)
biomet$timestamp <- ts

# If timestamp still has NA's - fill them in manually according to the closest previous half-hour
datasets <- list(eddypro = eddypro, biomet = biomet)

for (nm in names(datasets)) {
  df <- datasets[[nm]]
  ts <- df$timestamp
  
  if (!any(is.na(ts))) { assign(nm, df, inherits=TRUE); next }
  
  for (i in seq_along(ts)) {
    if (is.na(ts[i]) && i > 1 && !is.na(ts[i-1])) {
      # step in UTC to avoid DST gaps, then view back in UTC
      prev_utc <- with_tz(ts[i-1], "UTC")
      fill_utc <- prev_utc + seconds(1800)
      ts[i]    <- with_tz(fill_utc, "UTC")
    }
  }
  
  df$timestamp <- ts
  datasets[[nm]] <- df
  assign(nm, df, inherits = TRUE)
}

#---------------------------------------------------------
# Some changes to eddypro and biomet

# Adding/renaming variables, e.g. MO_LENGTH and V_SIGMA
eddypro <- eddypro %>%
  rename(MO_LENGTH_fluxnet = L) %>%
  mutate(V_SIGMA_fluxnet = sqrt(v_var), air_temperature  = ifelse(air_temperature >= -30 + 273.15 & air_temperature <= 40 + 273.15, air_temperature, NA_real_))

# Extra columns for meta
extra <- read_excel(extra_path)

extra_ts <- extra %>%
  mutate(timestamp = as.POSIXct(paste(date, "23:30:00"), tz = "UTC")) %>%
  select(-date)

hh_seq <- data.frame(
  timestamp = seq(
    from = as.POSIXct("2020-01-01 00:00:00", tz = "UTC"),
    to   = as.POSIXct("2020-12-31 23:30:00", tz = "UTC"),
    by   = "30 min"
  )
)

extra <- hh_seq %>%
  left_join(extra_ts, by = "timestamp") %>%
  fill(canopy_height, displacement_height, roughness_height, .direction = "downup")

eddypro <- hh_seq %>%
  left_join(eddypro, by = "timestamp") %>%
  left_join(extra, by = "timestamp") %>%
  rename(displacement_height_meta = displacement_height, roughness_length_meta = roughness_height) %>%
  mutate(master_sonic_height_meta = 2.2, latitude_meta = 52.2982) 

# Biomet changes
biomet <- biomet %>%
  rename(PPFD = `PPFD (µmol m-2 s-1)`, AirTemp = `Air Temperature (°C)`, soil_water_T_avg_1 = `soil_water_T_Avg(1)`, soil_water_T_avg_2 = `soil_water_T_Avg(2)`) %>%
  mutate(Rg = PPFD / 2.0565)

# Change air_temperature
temp <- biomet %>%
  dplyr::select(timestamp, AirTemp) %>%
  rename(air_temperature = AirTemp)

eddypro <- eddypro %>%
  dplyr::select(-air_temperature) %>%
  left_join(temp, by = "timestamp")

#---------------------------------------------------------
# Deduplicate timestamps
eddypro <- eddypro %>% arrange(timestamp) %>% distinct(timestamp, .keep_all = TRUE)
biomet  <- biomet  %>% arrange(timestamp) %>% distinct(timestamp, .keep_all = TRUE)

#----------------------------------------------------------
# Merge
eddypro  <- eddypro  %>% filter(year(timestamp) == year_no)
biomet  <- biomet  %>% filter(year(timestamp) == year_no) %>% rename_with(~ paste0(.x, "_biomet"),  -timestamp)

merged <- eddypro %>%
  left_join(biomet,  by = "timestamp") %>%
  mutate(across(where(is.numeric), ~ na_if(.x, -9999))) %>%
  rename(u_star = `u*`) %>%
  mutate(air_temperature = na.approx(air_temperature, x = timestamp, rule =2), PPFD_biomet = na.approx(PPFD_biomet, x = timestamp, rule =2))

# Add AGC filter and disallow night-time uptake
agc <- read_excel(agc_path)
merged <- merged %>%
  left_join(agc, by = "timestamp") %>%
  mutate(co2_flux = ifelse(agc < th_agc, co2_flux, NA)) %>%
  mutate(co2_flux = ifelse(co2_flux < 0 & PPFD_biomet < 10, NA, co2_flux))

#----------------------------------------------------------
# Base filters 
merged <- merged %>%
  mutate(
    # 1) CO₂ range
    co2_flux_low_high = ifelse(
      replace_na(co2_flux < th_co2_low | co2_flux > th_co2_high, FALSE),
      NA_real_, as.numeric(co2_flux)
    ),
    
    # # 2) Signal strength (use only this column)
    # co2_flux_sig_strength = ifelse(
    #   replace_na(co2_signal_strength_7200_mean < th_sig_strength, FALSE),
    #   NA_real_, as.numeric(co2_flux)
    # ),
    
    # 3) Random error (use rand_err_co2_flux)
    co2_flux_rand_err = ifelse(
      replace_na(as.numeric(rand_err_co2_flux) > th_rand_err, FALSE),
      NA_real_, as.numeric(co2_flux)
    )
    
    # # 4) Flow rate range (CUSTOM_FLOWRATE_MEAN_fluxnet)
    # co2_flux_flow_range = ifelse(
    #   replace_na(
    #     (as.numeric(BADM_INST_GA_CP_TUBE_FLOW_RATE_GA_CO2_fluxnet) < th_flow_low) |
    #       (as.numeric(BADM_INST_GA_CP_TUBE_FLOW_RATE_GA_CO2_fluxnet) > th_flow_high),
    #     FALSE
    #   ),
    #   NA_real_, as.numeric(co2_flux)
    # ) 
    
    # # 5) PPFD low-light masking ONLY when flux is negative
    # ,co2_flux_ppfd = ifelse(
    #   replace_na( (as.numeric(PPFD_1_1_1_biomet) < th_ppfd_day) & (as.numeric(co2_flux) < 0), FALSE ),
    #   NA_real_, as.numeric(co2_flux)
    # ),
    # 
    # # 6) u* filter
    # co2_flux_ustar = ifelse(
    #   replace_na(as.numeric(u_star) < th_ustar, FALSE),
    #   NA_real_, as.numeric(co2_flux)
    # )
  )

# which filtered columns to plot (adjust if you added/renamed)
filter_cols <- c(
  "co2_flux_low_high",
  # "co2_flux_sig_strength",
  "co2_flux_rand_err"
  # "co2_flux_flow_range"
  #, "co2_flux_ppfd",
  # "co2_flux_ustar"
)

plot_df <- merged %>%
  pivot_longer(
    cols = any_of(filter_cols),
    names_to = "filter",
    values_to = "co2_flux_filtered"
  ) %>%
  filter(!is.na(co2_flux_filtered)) %>%
  mutate(
    filter = recode(filter,
                    co2_flux_low_high      = "Low/High range",
                    # co2_flux_sig_strength  = "Signal strength",
                    co2_flux_rand_err      = "Random error"
                    # co2_flux_flow_range    = "Flow range"
                    #, co2_flux_ppfd          = "PPFD (low & flux<0)",
                    #co2_flux_ustar         = "u*"
    )
  )

ggplot(plot_df, aes(x = timestamp, y = co2_flux_filtered)) +
  geom_point(size = 0.6, alpha = 0.6, na.rm = TRUE) +
  scale_x_datetime(
    date_breaks = "7 days",
    date_labels = "%Y/%m/%d",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(x = NULL, y = "CO₂ flux", title = "CO₂ flux with single filters applied") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap(~ filter, ncol = 2, scales = "free_y")

# Apply all base filters
merged <- merged %>%
  mutate(
    fail_base =
      replace_na(co2_flux < th_co2_low  | co2_flux > th_co2_high, FALSE) |
      # replace_na(co2_signal_strength_7200_mean < th_sig_strength, FALSE) |
      replace_na(as.numeric(rand_err_co2_flux) > th_rand_err, FALSE),
    # replace_na( (as.numeric(PPFD_1_1_1_biomet) < th_ppfd_day) & (as.numeric(co2_flux) < 0), FALSE ) |
    # replace_na(as.numeric(u_star) < th_ustar, FALSE),
    
    co2_flux_base_filters = ifelse(fail_base, NA_real_, as.numeric(co2_flux))
  ) %>%
  select(-fail_base)

#----------------------------------------------------------
# Seasonal u* threshold (REddyProc)
ep_df <- merged %>%
  transmute(
    DateTime = timestamp,   
    Year  = year(timestamp),
    DoY   = yday(timestamp),
    Hour  = hour(timestamp) + minute(timestamp) / 60,
    NEE   = ifelse(is.finite(co2_flux_base_filters),        co2_flux_base_filters,        NA_real_),
    Ustar = ifelse(is.finite(u_star),           u_star,          NA_real_),
    Tair  = ifelse(is.finite(air_temperature),
                   air_temperature - 273.15,                     # Kelvin → Celsius
                   NA_real_),
    Rg    = ifelse(is.finite(Rg_biomet),               Rg_biomet,              NA_real_)
  )

hh_skeleton <- data.frame(
  DateTime = seq(
    from = as.POSIXct("2020-01-01 00:30:00", tz = "UTC"),
    to   = as.POSIXct("2020-12-31 23:30:00", tz = "UTC"),
    by   = "30 min"
  )
)

ep_df <- hh_skeleton %>%
  left_join(ep_df, by = "DateTime") %>%
  mutate(
    Year = year(DateTime),
    DoY  = yday(DateTime),
    Hour = hour(DateTime) + minute(DateTime) / 60
  )

season_factor <- usCreateSeasonFactorMonth(
  ep_df$DateTime,
  startMonth = c(2, 5, 8, 11)
)

# count actual valid nighttime records per season to set thresholds
n_per_season <- ep_df %>%
  mutate(season = season_factor) %>%
  filter(is.finite(NEE), is.finite(Tair), is.finite(Ustar), is.finite(Rg), Rg < 10) %>%
  count(season)

n_year   <- sum(n_per_season$n)
min_s    <- max(floor(min(n_per_season$n) * 0.5), 10)
min_t    <- max(floor(min_s / 7), 5)
min_y    <- max(n_year, 50)

ctrl_sub <- usControlUstarSubsetting(
  minRecordsWithinSeason = min_s,
  minRecordsWithinTemp   = min_t,
  minRecordsWithinYear   = min_y
)

EProc <- sEddyProc$new('JC1_2020', ep_df, c('NEE', 'Ustar', 'Tair', 'Rg'))
EProc$sEstUstarThold(seasonFactor = season_factor, ctrlUstarSub = ctrl_sub)

ustar_th <- EProc$sUSTAR_DETAILS$uStarTh %>%
  filter(aggregationMode == "season") %>%
  select(season, uStar_thresh = uStar)

season_factor_merged <- usCreateSeasonFactorMonth(
  merged$timestamp,
  startMonth = c(2, 5, 8, 11)
)

merged <- merged %>%
  mutate(
    ustar_season          = as.character(season_factor_merged),
    ustar_seasonal_thresh = ustar_th$uStar_thresh[
      match(ustar_season, as.character(ustar_th$season))
    ],
    co2_flux_base_filters = ifelse(
      replace_na(u_star < ustar_seasonal_thresh, FALSE),
      NA_real_, co2_flux_base_filters
    )
  )

# plot co2_flux_base_filters
ggplot(merged, aes(x = timestamp, y = co2_flux_base_filters)) +
  geom_point(size = 0.6, alpha = 0.6, na.rm = TRUE) +
  scale_x_datetime(
    date_breaks = "7 days",
    date_labels = "%Y/%m/%d",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(x = NULL, y = "CO₂ flux (base filters)", title = "CO₂ flux after base filters") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )
#----------------------------------------------------------
# QC flag filters (1-9)

for (i in 1:9) {
  nm <- paste0("co2_flux_base_filters_", i)
  merged <- merged %>%
    mutate(!!nm := ifelse(qc_co2_flux > i, NA_real_, co2_flux_base_filters))
}

# Plot qced fluxes
for (i in 1:9) {
  col <- paste0("co2_flux_base_filters_", i)
  
  p <- ggplot(merged, aes(x = timestamp, y = .data[[col]])) +
    geom_point(size = 0.6, alpha = 0.6, na.rm = TRUE) +
    scale_x_datetime(
      date_breaks = "7 days",
      date_labels = "%Y/%m/%d",
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = NULL,
      y = "CO₂ flux (base filters)",
      title = paste0("CO₂ flux after base filters — i = ", i)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  # ggsave(sprintf("co2_flux_base_filters_%d.png", i), p, width = 10, height = 5, dpi = 150)
}

#----------------------------------------------------------
# Grid footprint filter (Kljun 2015 2-D FFP inside the site polygon)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
  library(tibble)
})

# source the PBLH and 2-D footprint engines
source(here::here("data_preparation/quality_control_nee/grid_qc/pblh_calc.R"))
source(here::here("data_preparation/quality_control_nee/grid_qc/calc_footprint_FFP_mod.R"))

#----------------------------------------------------------
# boundary polygon + tower location (pick the right site row)
boundary_polygons_path <- here::here("data_preparation/quality_control_nee/grid_qc/nasco_site_list_wkt.csv")
poly_tbl <- readr::read_csv(boundary_polygons_path, show_col_types = FALSE)

row_poly <- poly_tbl %>%
  filter(`Site Name` == site_key) %>%
  slice(1)

boundary_poly <- sf::st_as_sfc(row_poly$boundary_2157, crs = 2157)
tower_pt      <- sf::st_as_sfc(row_poly$tower_loc_2157, crs = 2157)
tower_xy      <- sf::st_coordinates(tower_pt)[1, ]
tower_x <- as.numeric(tower_xy[1])
tower_y <- as.numeric(tower_xy[2])

#----------------------------------------------------------
# inputs required for footprint (create/repair when missing)

# measurement height above displacement
if (!"zm" %in% names(merged)) {
  merged <- merged %>%
    mutate(zm = master_sonic_height_meta - displacement_height_meta)
}

# fallback air_temperature column if absent to avoid errors
if (!"air_temperature" %in% names(merged)) {
  merged$air_temperature <- NA_real_
}

# (z-d)/L needed by your pblh function; if missing, approximate as zm / MO_LENGTH_fluxnet
if (!("(z-d)/L" %in% names(merged))) {
  merged <- merged %>%
    mutate(`(z-d)/L` = zm / MO_LENGTH_fluxnet)
}

# latitude fallback (52 deg) if missing
merged <- merged %>%
  mutate(lat_for_pblh = dplyr::coalesce(latitude_meta, 52))

#----------------------------------------------------------
# pblh computed row-by-row (avoid vectorised internal state issues)

pblh_vec <- rep(NA_real_, nrow(merged))

ok_pblh <- with(merged,
                is.finite(MO_LENGTH_fluxnet) &
                  is.finite(u_star) &
                  is.finite(V_SIGMA_fluxnet) &
                  is.finite(lat_for_pblh) &
                  is.finite(`(z-d)/L`) &
                  is.finite(air_temperature)
)

idx <- which(ok_pblh)

for (i in idx) {
  pblh_vec[i] <- tryCatch(
    boundary_layer_height(
      Ls     = merged$MO_LENGTH_fluxnet[i],
      ustars = merged$u_star[i],
      t_covs = merged$V_SIGMA_fluxnet[i],
      LAT    = merged$lat_for_pblh[i],
      zLs    = merged$`(z-d)/L`[i],
      air_ts = merged$air_temperature[i]
    )[1],
    error = function(e) NA_real_
  )
}

merged$pblh <- pblh_vec

#----------------------------------------------------------
# compute footprint ratio per row (flux fraction inside boundary polygon)

req_cols <- c("zm","roughness_length_meta","wind_speed","pblh",
              "MO_LENGTH_fluxnet","V_SIGMA_fluxnet","u_star","wind_dir")

ok_fp <- merged %>%
  transmute(across(all_of(req_cols), ~ is.finite(.))) %>%
  mutate(all_ok = rowSums(across(everything())) == length(req_cols)) %>%
  pull(all_ok)

n <- nrow(merged)
footprint_ratio <- rep(NA_real_, n)

rows_to_run <- which(ok_fp)

for (i in rows_to_run) {
  fp <- try({
    calc_footprint_FFP_mod(
      zm       = merged$zm[i],
      z0       = merged$roughness_length_meta[i],
      umean    = merged$wind_speed[i],
      h        = merged$pblh[i],
      ol       = merged$MO_LENGTH_fluxnet[i],
      sigmav   = merged$V_SIGMA_fluxnet[i],
      ustar    = merged$u_star[i],
      wind_dir = merged$wind_dir[i],
      nx       = 200,
      tower_x  = tower_x,
      tower_y  = tower_y
    )
  }, silent = TRUE)
  
  if (inherits(fp, "try-error") || is.null(fp$f_2d)) next
  
  # fraction of footprint mass within boundary polygon
  x_vec <- as.vector(fp$x_2d)
  y_vec <- as.vector(fp$y_2d)
  f_vec <- as.vector(fp$f_2d)
  
  pts <- sf::st_as_sf(data.frame(x = x_vec, y = y_vec), coords = c("x","y"), crs = 2157)
  inside <- sf::st_within(pts, boundary_poly, sparse = FALSE)[,1]
  
  tot <- sum(f_vec, na.rm = TRUE)
  inb <- sum(f_vec[inside], na.rm = TRUE)
  
  footprint_ratio[i] <- ifelse(is.finite(tot) && tot > 0, inb / tot, NA_real_)
}

merged$footprint_ratio_grid <- footprint_ratio

# binary footprint qc flags
merged$flux_qc_footprint_grid_70 <- as.integer(merged$footprint_ratio_grid >= 0.70)
merged$flux_qc_footprint_grid_80 <- as.integer(merged$footprint_ratio_grid >= 0.80)
merged$flux_qc_footprint_grid_90 <- as.integer(merged$footprint_ratio_grid >= 0.90)

#----------------------------------------------------------
# apply grid footprint qc to each co2_flux_base_filters_i (i = 1..9)

for (i in 1:9) {
  src   <- paste0("co2_flux_base_filters_", i)
  
  out70 <- paste0("co2_flux_base_filters_", i, "_70_grid")
  out80 <- paste0("co2_flux_base_filters_", i, "_80_grid")
  out90 <- paste0("co2_flux_base_filters_", i, "_90_grid")
  
  merged[[out70]] <- ifelse(merged$flux_qc_footprint_grid_70 == 1L, merged[[src]], NA_real_)
  merged[[out80]] <- ifelse(merged$flux_qc_footprint_grid_80 == 1L, merged[[src]], NA_real_)
  merged[[out90]] <- ifelse(merged$flux_qc_footprint_grid_90 == 1L, merged[[src]], NA_real_)
}

# save as CSV
write_csv(merged, file = file.path(results_path, "merged_qc.csv"))

# or save as RDS (keeps types exactly, faster to reload)
saveRDS(merged, file = file.path(results_path, "merged_qc.rds"))

