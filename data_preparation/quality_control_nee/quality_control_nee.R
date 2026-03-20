# =============================================================================
# quality_control_nee.R  — NEE Quality Control Pipeline
# =============================================================================
#
# PURPOSE
#   Applies a multi-stage quality control (QC) procedure to raw CO2 flux data
#   from EddyPro, producing a merged, flagged, and footprint-filtered dataset
#   that can be used directly by the data preparation pipeline (Script 02+).
#
# DATA SOURCES
#   eddypro  -- EddyPro full-output CSV: co2_flux, qc_co2_flux, u_star,
#               wind_dir, wind_speed, MO_LENGTH, rand_err_co2_flux, V_SIGMA,
#               co2_signal_strength_7200_mean, date, time
#   biomet   -- Biomet half-hourly CSV: PPFD_1_1_1, sonic height, etc.
#   meta     -- Site metadata CSV: sonic height, displacement height,
#               roughness length, latitude
#   fluxnet  -- FluxNet-format CSV: TIMESTAMP_START, flow rate, MO_LENGTH,
#               FETCH_70/80/90, V_SIGMA
#
# PIPELINE STAGES
#   1. Paths and thresholds  -- set input file paths and filter thresholds
#   2. Load data             -- read CSV or XLSX files auto-detected by extension
#   3. Timestamps            -- parse and repair all four datasets
#   4. Merge                 -- left-join all sources onto the EddyPro base;
#                               replace -9999 sentinels with NA
#   5. Base filters          -- apply 4 independent QC criteria, each producing
#                               a separate filtered flux column for inspection
#   6. Combined base filter  -- apply all criteria simultaneously into
#                               co2_flux_base_filters
#   7. Stationarity flag     -- apply EddyPro qc_co2_flux flags 1-9, producing
#                               co2_flux_base_filters_{1..9}
#   8. Grid footprint QC     -- compute PBLH and Kljun FFP footprint fraction;
#                               apply 70/80/90% boundary thresholds, producing
#                               co2_flux_base_filters_{1..9}_{70|80|90}_grid
#   9. Save outputs          -- write merged_qc.csv and merged_qc.rds
#
# OUTPUT COLUMN NAMING CONVENTION
#   co2_flux_base_filters         -- base filters only (no stationarity flag)
#   co2_flux_base_filters_{i}     -- base + stationarity flag <= i (i = 1..9)
#   co2_flux_base_filters_{i}_70_grid -- above + footprint >= 70% inside boundary
#   co2_flux_base_filters_{i}_80_grid -- above + footprint >= 80%
#   co2_flux_base_filters_{i}_90_grid -- above + footprint >= 90%
#   (Script 03 in the data preparation pipeline reads the _6_70_grid column.)
#
# OUTPUTS  (written to results_path/)
#   merged_qc.csv   -- complete merged and flagged data frame
#   merged_qc.rds   -- same data in RDS format (faster to reload, preserves types)
#
# HOW TO RUN
#   1. Fill in the four *_path variables in Section 1.
#   2. Set results_path to the desired output directory.
#   3. Review and adjust the filter thresholds in Section 1.
#   4. Source this file:  source("quality_control_nee.R")
#   The footprint section (Section 8) requires the site boundary polygon CSV
#   and tower location; update site_key to match your site.
#
# =============================================================================

# =============================================================================
# SECTION 1 — Packages
# =============================================================================
library(here)
library(readr)
library(readxl)
library(tidyverse)
library(hms)
# =============================================================================
# SECTION 2 — Paths, output directory, and filter thresholds
# =============================================================================
# --- File paths (fill in before running) ------------------------------------
# Each path accepts either a .csv or an .xlsx/.xls file; the correct reader
# is chosen automatically in Section 3.
eddypro_path <- here::here("")
meta_path    <- here::here("")
fluxnet_path <- here::here("")
biomet_path  <- here::here("")

# Output directory -- rename to results_qc_{SITE}_{YEAR} before running
results_path <- here::here("results_qc_site_year")
dir.create(results_path, recursive = TRUE, showWarnings = FALSE)

# --- Filter thresholds -------------------------------------------------------
# th_co2_low / th_co2_high: physiologically plausible CO2 flux range.
#   Values outside this window are almost always instrument artefacts.
#   Adjust for very productive or very disturbed ecosystems if needed.
# th_rand_err: maximum acceptable random error (umol m-2 s-1).
#   Large random errors indicate unstable flux calculation conditions.
# th_sig_strength: minimum acceptable LI-7200 CO2 signal strength (%).
#   Values below this indicate a dirty or damp optical cell.
# th_flow_low / th_flow_high: valid LI-7200 sample tube flow rate (L/min).
#   Outside-range flow indicates pump failure or a blocked inlet tube.
# Commented thresholds (th_ppfd_day, th_ustar) are available if needed;
#   uncomment and activate the corresponding filter blocks below.
th_co2_low      <- -40
th_co2_high     <-  30
th_fetch70      <- NA_real_   # e.g. 100 to require FETCH_70_fluxnet >= 100
th_fetch80      <- NA_real_
th_fetch90      <- NA_real_
th_rand_err     <- 100        # max random error
th_sig_strength <- 70         # min signal strength
th_flow_low     <- 12         # valid flow range low
th_flow_high    <- 18         # valid flow range high
# th_ppfd_day     <- 10         # ppfd daytime threshold
# th_ustar        <- 0.1

# =============================================================================
# SECTION 3 — Load data
# =============================================================================
# Each source file is read with either readr::read_csv() or readxl::read_excel()
# depending on the file extension.  The extension check is case-insensitive.
# An informative error is raised for any unsupported format.

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

## meta
if (grepl("\\.csv$", meta_path,  ignore.case = TRUE)) {
  meta <- readr::read_csv(meta_path)
} else if (grepl("\\.xlsx?$", meta_path, ignore.case = TRUE)) {
  meta <- readxl::read_excel(meta_path)
} else stop("Unsupported file type for meta_path")

## fluxnet
if (grepl("\\.csv$", fluxnet_path,  ignore.case = TRUE)) {
  fluxnet <- readr::read_csv(fluxnet_path)
} else if (grepl("\\.xlsx?$", fluxnet_path, ignore.case = TRUE)) {
  fluxnet <- readxl::read_excel(fluxnet_path)
} else stop("Unsupported file type for fluxnet_path")

# =============================================================================
# SECTION 4 — Timestamp construction and repair
# =============================================================================
# EddyPro / biomet / meta: timestamps are built from two separate columns
# ('date' and 'time') that are present in the raw output.  lubridate::
# ymd_hms() parses the 'time' string and update() replaces its date
# components with those from the 'date' column.  Timezone is set to
# Europe/Dublin so that DST transitions are handled correctly; all four
# files are then harmonised to the same tz before the merge.
#
# FluxNet: timestamp is encoded as an integer in YYYYMMDDHHMM format;
# ymd_hm(as.numeric(...)) parses it directly.
#
# NA repair: if any timestamp is still NA after parsing (e.g. due to a
# malformed row in the raw file), it is filled by stepping forward 30 min
# from the most recent non-NA predecessor.  The step is performed in UTC
# to avoid DST jumps, then converted back to Europe/Dublin.
# A final anyDuplicated() check ensures no two rows share the same timestamp
# before the merge.

eddypro <- eddypro %>%
  mutate(
    timestamp = update(ymd_hms(time, tz = "Europe/Dublin"),
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
    timestamp = update(ymd_hms(time, tz = "Europe/Dublin"),
                       year  = year(as_date(date)),
                       month = month(as_date(date)),
                       mday  = mday(as_date(date)))
  )  %>%
  relocate(timestamp, .before = 1)

ts <- biomet$timestamp
for (i in seq_along(ts)) if (i > 1 && is.na(ts[i]) && !is.na(ts[i-1])) ts[i] <- ts[i-1] + minutes(30)
biomet$timestamp <- ts

meta <- meta %>%
  mutate(
    timestamp = update(ymd_hms(time, tz = "Europe/Dublin"),
                       year  = year(as_date(date)),
                       month = month(as_date(date)),
                       mday  = mday(as_date(date)))
  )  %>%
  relocate(timestamp, .before = 1)

ts <- meta$timestamp
for (i in seq_along(ts)) if (i > 1 && is.na(ts[i]) && !is.na(ts[i-1])) ts[i] <- ts[i-1] + minutes(30)
meta$timestamp <- ts

fluxnet <- fluxnet %>%
  mutate(timestamp = ymd_hm(as.numeric(TIMESTAMP_START)))  %>%
  relocate(timestamp, .before = 1)

# If timestamp still has NA's - fill them in manually according to the closest previous half-hour
datasets <- list(eddypro = eddypro, biomet = biomet, meta = meta, fluxnet = fluxnet)

for (nm in names(datasets)) {
  df <- datasets[[nm]]
  ts <- df$timestamp
  
  if (!any(is.na(ts))) { cat(nm, ": no NA timestamps\n", sep=""); assign(nm, df, inherits=TRUE); next }
  
  cat("\n", nm, " — filling NA timestamps (30-min steps via UTC):\n", sep="")
  for (i in seq_along(ts)) {
    if (is.na(ts[i]) && i > 1 && !is.na(ts[i-1])) {
      # step in UTC to avoid DST gaps, then view back in Europe/Dublin
      prev_utc <- with_tz(ts[i-1], "UTC")
      fill_utc <- prev_utc + seconds(1800)
      ts[i]    <- with_tz(fill_utc, "Europe/Dublin")
      cat(sprintf("  row %d: set to %s (from %s + 30 min)\n",
                  i,
                  format(ts[i], "%Y-%m-%d %H:%M:%S %Z"),
                  format(ts[i-1], "%Y-%m-%d %H:%M:%S %Z")))
    }
  }
  
  df$timestamp <- ts
  datasets[[nm]] <- df
  assign(nm, df, inherits = TRUE)
}

# Check remaining NAs
sapply(datasets, \(x) sum(is.na(x$timestamp)))



# You will need to fill timestamps yourself if stll missing
na_timestamp <- sapply(list(eddypro=eddypro, biomet=biomet, meta=meta, fluxnet = fluxnet), \(x) sum(is.na(x$timestamp)))
print(na_timestamp)

if (any(na_timestamp > 0)) {
  message("Some NA values remain in timestamp — please fill them in manually.")
} else {
  message("ll timestamp values are filled — no manual edits needed.")
}

# Deduplicate timestamps
eddypro <- eddypro %>% arrange(timestamp) %>% distinct(timestamp, .keep_all = TRUE)
biomet  <- biomet  %>% arrange(timestamp) %>% distinct(timestamp, .keep_all = TRUE)
meta    <- meta    %>% arrange(timestamp) %>% distinct(timestamp, .keep_all = TRUE)
fluxnet <- fluxnet    %>% arrange(timestamp) %>% distinct(timestamp, .keep_all = TRUE)

# Check for duplicates in timestamp
anyDuplicated(eddypro$timestamp)
anyDuplicated(biomet$timestamp)
anyDuplicated(meta$timestamp)
anyDuplicated(fluxnet$timestamp)

# =============================================================================
# SECTION 5 — Merge all sources onto EddyPro base
# =============================================================================
# Before merging, each secondary source is renamed with a suffix (_biomet,
# _fluxnet, _meta) to prevent column-name collisions when all four data frames
# are joined.  EddyPro columns retain their original names as the primary source.
# The join uses left_join on 'timestamp', so all EddyPro rows are retained;
# biomet/meta/fluxnet rows without a matching EddyPro timestamp are dropped.
# The -9999 EddyPro sentinel (missing value indicator) is replaced with NA
# across all numeric columns so that downstream filters see a consistent NA.
# u* is renamed to u_star to avoid backtick syntax throughout the rest of the script.

biomet  <- biomet  %>% rename_with(~ paste0(.x, "_biomet"),  -timestamp)
fluxnet <- fluxnet %>% rename_with(~ paste0(.x, "_fluxnet"), -timestamp)
meta    <- meta    %>% rename_with(~ paste0(.x, "_meta"),    -timestamp)

merged <- eddypro %>%
  left_join(biomet,  by = "timestamp") %>%
  left_join(meta,    by = "timestamp") %>%
  left_join(fluxnet, by = "timestamp") %>%
  mutate(across(where(is.numeric), ~ na_if(.x, -9999))) %>%
  rename(u_star = `u*`)

# =============================================================================
# SECTION 6 — Individual base filters (diagnostic)
# =============================================================================
# Four independent QC criteria are applied, each producing a separate filtered
# flux column.  These columns are plotted individually so that the operator can
# judge how much data each filter removes before committing to the combined filter
# in Section 7.  Filters 5 (PPFD masking) and 6 (u* threshold) are commented out;
# uncomment them if low-turbulence or low-light exclusion is required.
#
#   Filter 1 -- CO2 range: removes physically implausible flux values.
#   Filter 2 -- Signal strength: removes half-hours with a dirty/wet optical cell.
#   Filter 3 -- Random error: removes half-hours with unstable flux estimates.
#   Filter 4 -- Flow rate: removes half-hours when the sample tube flow is abnormal.
merged <- merged %>%
  mutate(
    # 1) CO2 range
    co2_flux_low_high = ifelse(
      replace_na(co2_flux < th_co2_low | co2_flux > th_co2_high, FALSE),
      NA_real_, as.numeric(co2_flux)
    ),
    
    # 2) Signal strength (use only this column)
    co2_flux_sig_strength = ifelse(
      replace_na(co2_signal_strength_7200_mean < th_sig_strength, FALSE),
      NA_real_, as.numeric(co2_flux)
    ),
    
    # 3) Random error (use rand_err_co2_flux)
    co2_flux_rand_err = ifelse(
      replace_na(as.numeric(rand_err_co2_flux) > th_rand_err, FALSE),
      NA_real_, as.numeric(co2_flux)
    ),
    
    # 4) Flow rate range (CUSTOM_FLOWRATE_MEAN_fluxnet)
    co2_flux_flow_range = ifelse(
      replace_na(
        (as.numeric(BADM_INST_GA_CP_TUBE_FLOW_RATE_GA_CO2_fluxnet) < th_flow_low) |
          (as.numeric(BADM_INST_GA_CP_TUBE_FLOW_RATE_GA_CO2_fluxnet) > th_flow_high),
        FALSE
      ),
      NA_real_, as.numeric(co2_flux)
    ) 
    
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
  "co2_flux_sig_strength",
  "co2_flux_rand_err",
  "co2_flux_flow_range"
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
                    co2_flux_sig_strength  = "Signal strength",
                    co2_flux_rand_err      = "Random error",
                    co2_flux_flow_range    = "Flow range"
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


# =============================================================================
# SECTION 7 — Combined base filter
# =============================================================================
# Applies all active base filter criteria simultaneously (logical OR of all
# fail conditions).  Any half-hour that fails at least one criterion is set
# to NA in co2_flux_base_filters.  This is the canonical 'passed base QC'
# column and the starting point for all subsequent stationarity and footprint
# filters.

merged <- merged %>%
  mutate(
    fail_base =
      replace_na(co2_flux < th_co2_low  | co2_flux > th_co2_high, FALSE) |
      replace_na(co2_signal_strength_7200_mean < th_sig_strength, FALSE) |
      replace_na(as.numeric(rand_err_co2_flux) > th_rand_err, FALSE) |
      replace_na(as.numeric(BADM_INST_GA_CP_TUBE_FLOW_RATE_GA_CO2_fluxnet) < th_flow_low,  FALSE) |
      replace_na(as.numeric(BADM_INST_GA_CP_TUBE_FLOW_RATE_GA_CO2_fluxnet) > th_flow_high, FALSE),
    # replace_na( (as.numeric(PPFD_1_1_1_biomet) < th_ppfd_day) & (as.numeric(co2_flux) < 0), FALSE ) |
    # replace_na(as.numeric(u_star) < th_ustar, FALSE),
    
    co2_flux_base_filters = ifelse(fail_base, NA_real_, as.numeric(co2_flux))
  ) %>%
  select(-fail_base)

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
# =============================================================================
# SECTION 8 — Stationarity and developed-turbulence flag filters (qc_co2_flux)
# =============================================================================
# EddyPro's qc_co2_flux column encodes a combined stationarity and integral
# turbulence characteristic test, reported on a 0-9 scale (Foken et al., 2004):
#   0-1  = best quality (fully developed turbulence, stationary)
#   2-6  = moderate quality
#   7-9  = poor quality
# For each threshold i, co2_flux_base_filters_{i} retains only half-hours where
# qc_co2_flux <= i.  The standard selection for long-term flux studies is i = 6
# (co2_flux_base_filters_6); the downstream pipeline reads the _6_70_grid column.
# Plotting all nine thresholds below allows the operator to choose a different
# cutoff if the distribution of quality flags at this site warrants it.

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

# =============================================================================
# SECTION 9 — Grid-based flux footprint quality control
# =============================================================================
# Computes the fraction of the flux footprint that falls inside the site
# boundary polygon at each half-hour, and applies threshold filters at 70%,
# 80%, and 90%.  Half-hours where the footprint is predominantly outside the
# managed grassland boundary are excluded to prevent measurements of adjacent
# land cover types (arable, hedgerow, road) from contaminating the NEE record.
#
# ALGORITHM
#   Step 1 -- Compute planetary boundary layer height (PBLH) per half-hour
#             using boundary_layer_height() (pblh_calc.R).
#             PBLH is required by the Kljun footprint model as an upper
#             boundary condition on the turbulent mixing layer.
#   Step 2 -- For each half-hour with complete inputs, call
#             calc_footprint_FFP_mod() (Kljun et al., 2015; Flux Footprint
#             Prediction, FFP) to generate a 2D probability density grid
#             in Irish Transverse Mercator coordinates (EPSG:2157).
#   Step 3 -- Count what fraction of the total footprint flux mass falls
#             inside the site boundary polygon using sf::st_within().
#   Step 4 -- Apply binary threshold flags (footprint >= 70/80/90%) and
#             mask the flux columns accordingly.
#
# INPUTS
#   grid_qc/nasco_site_list_wkt.csv  -- site boundary polygons in WKT (EPSG:2157)
#   grid_qc/pblh_calc.R              -- boundary_layer_height() function
#   grid_qc/calc_footprint_FFP_mod.R -- calc_footprint_FFP_mod() function
#
# REFERENCE
#   Kljun, N., Calanca, P., Rotach, M.W., Schmid, H.P. (2015). A simple
#   two-dimensional parameterisation for Flux Footprint Prediction (FFP).
#   Geoscientific Model Development, 8, 3695-3713.


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(sf)
  library(tibble)
})

# kljun + pblh engines (keep your originals)
source(here::here("grid_qc/pblh_calc.R"))             # boundary_layer_height()
source(here::here("grid_qc/calc_footprint_FFP_mod.R"))# calc_footprint_FFP_mod()

#----------------------------------------------------------
# boundary polygon + tower location (pick the right site row)
boundary_polygons_path <- here::here("grid_qc/nasco_site_list_wkt.csv")
poly_tbl <- readr::read_csv(boundary_polygons_path, show_col_types = FALSE)

# Change site_key to the 'Site Name' value in nasco_site_list_wkt.csv that
# corresponds to the current site.  The boundary polygon and tower coordinates
# are extracted from this row.
site_key <- "Johnstown Castle - Teagasc (1)"   # JC1; use (2) for JC2

row_poly <- poly_tbl %>%
  filter(`Site Name` == site_key) %>%
  slice(1)

boundary_poly <- sf::st_as_sfc(row_poly$boundary_2157, crs = 2157)
tower_pt      <- sf::st_as_sfc(row_poly$tower_loc_2157, crs = 2157)
tower_xy      <- sf::st_coordinates(tower_pt)[1, ]
tower_x <- as.numeric(tower_xy[1])
tower_y <- as.numeric(tower_xy[2])

# --- Required input columns (create or repair if absent in merged) -----------
# zm: effective measurement height = sonic height - displacement height.
#   Computed from meta columns; displacement height is typically ~2/3 of
#   canopy height and can vary seasonally in tall vegetation.
# air_temperature: used by PBLH unstable-conditions formula; set to NA_real_
#   if absent so rows with missing temperature get NA PBLH (not an error).
# (z-d)/L: Monin-Obukhov stability parameter; if not precomputed by EddyPro,
#   approximated here as zm / MO_LENGTH_fluxnet.
# lat_for_pblh: site latitude used for the Coriolis parameter in PBLH.
#   Defaults to 52 deg (central Ireland) if latitude_meta is absent.

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

# --- Compute PBLH row-by-row -------------------------------------------------
# boundary_layer_height() maintains an internal state variable (h_cur) that
# tracks the evolving boundary layer depth under unstable conditions.  This
# state depends on the order of rows and cannot be vectorised without
# replicating the sequential logic.  The loop is therefore kept explicit.
# Rows with any missing input receive NA PBLH (tryCatch catches edge cases
# such as a zero Coriolis parameter at very low latitudes).
# ok_pblh pre-filters rows with at least one non-finite input to avoid
# calling the function when it would trivially return NA.

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

# --- Compute per-half-hour footprint fraction --------------------------------
# For each row with complete inputs, calc_footprint_FFP_mod() returns a 2D
# probability density matrix (f_2d) on a regular x-y grid centred at the tower
# and rotated to the prevailing wind direction.  Grid coordinates are in the
# same projected CRS (EPSG:2157) as the boundary polygon.
# sf::st_within() identifies which grid cells fall inside the boundary; the
# footprint_ratio is the sum of f_2d values inside the boundary divided by
# the total sum (i.e. the footprint mass fraction inside the site).
# nx = 200 controls the spatial resolution of the footprint grid; increase
# for higher precision at the cost of computation time.

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

# --- Apply footprint threshold flags to all stationarity-filtered columns ----
# For each of the 9 stationarity-flag columns, three new columns are created
# by masking flux values where the footprint ratio is below 70%, 80%, or 90%.
# This produces a 3x9 = 27-column matrix of flux series at different combined
# QC stringency levels.  Script 03 reads co2_flux_base_filters_6_70_grid as
# the standard analysis column (stationarity <= 6, footprint >= 70%).

for (i in 1:9) {
  src   <- paste0("co2_flux_base_filters_", i)
  
  out70 <- paste0("co2_flux_base_filters_", i, "_70_grid")
  out80 <- paste0("co2_flux_base_filters_", i, "_80_grid")
  out90 <- paste0("co2_flux_base_filters_", i, "_90_grid")
  
  merged[[out70]] <- ifelse(merged$flux_qc_footprint_grid_70 == 1L, merged[[src]], NA_real_)
  merged[[out80]] <- ifelse(merged$flux_qc_footprint_grid_80 == 1L, merged[[src]], NA_real_)
  merged[[out90]] <- ifelse(merged$flux_qc_footprint_grid_90 == 1L, merged[[src]], NA_real_)
}

#----------------------------------------------------------
# quick retention glance (example uses i = 3; change i if you like)

i <- 3
ret_grid <- tibble::tibble(
  perc = c(70, 80, 90),
  kept = c(
    mean(!is.na(merged[[paste0("co2_flux_base_filters_", i, "_70_grid")]]), na.rm = TRUE),
    mean(!is.na(merged[[paste0("co2_flux_base_filters_", i, "_80_grid")]]), na.rm = TRUE),
    mean(!is.na(merged[[paste0("co2_flux_base_filters_", i, "_90_grid")]]), na.rm = TRUE)
  )
)
print(ret_grid)


# =============================================================================
# SECTION 10 — Save outputs
# =============================================================================
# Two output formats are saved:
#   merged_qc.csv  -- plain-text CSV for inspection in spreadsheet software
#                     and for use by Script 02 of the data preparation pipeline.
#   merged_qc.rds  -- R binary format; preserves all column types exactly
#                     and loads ~10x faster than CSV for large files.
# Script 02 reads merged_qc.csv (not the RDS) to allow non-R inspection.

write_csv(merged, file = file.path(results_path, "merged_qc.csv"))
saveRDS(merged,   file = file.path(results_path, "merged_qc.rds"))

