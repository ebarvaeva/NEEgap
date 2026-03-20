# =============================================================================
# ManagementEvents.R — Site Management Event Tables (SKELETON)
# =============================================================================
#
# PURPOSE
#   Defines all grazing and fertiliser event windows for every site and year,
#   then expands them into a standardised half-hourly event table that is used
#   by generate_management_csvs.R, 04_management_days_since.R, and
#   05_nitrogen.R.
#
# HOW TO FILL IN YOUR DATA
#   1. Add your grazing and fertiliser event dates to the `management` list
#      (Section 1) and to the per-site detail tables (Section 2).
#   2. The code_map (Section 0) translates farm-record grazing intensity
#      codes (0.2–1.0) into half-hour time windows.  Edit it only if your
#      recording system uses different window conventions.
#   3. Do NOT edit Section 3 (helper functions) or Section 4 (event building
#      and finalisation) — those are data-independent logic blocks.
#
# PLACEHOLDER CONVENTION
#   Every line that needs real data is marked:
#     "YYYY-MM-DD"    — replace with an actual date string
#     "DD.MM.YYYY"    — replace with a date in day.month.year format
#     0.0             — replace with the actual grazing intensity code (0.2–1.0)
#     0L              — replace with the actual animal count (integer)
#     0.0             — replace with the actual N amount (kg N ha⁻¹)
#
# SOURCED BY
#   generate_management_csvs.R (produces CSV files for the pipeline)
#   Regrowth_period.R          (adds regrowth period dummies to the RDS)
#
# =============================================================================

# Required packages
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(readr)
library(lubridate)
library(purrr)


# =============================================================================
# SECTION 0 — Code map
# =============================================================================
# Translates the farm-record grazing intensity code (0.2–1.0) into one, two,
# or three half-hour time windows for that day.  carry_next = TRUE means the
# grazing window continues into the first hours of the following day.
#
# Edit only if your field-recording system uses different time window
# conventions than those listed here.

code_map <- tibble::tribble(
  ~code, ~t1_start, ~t1_end, ~t2_start, ~t2_end, ~t3_start, ~t3_end, ~carry_next,
  1,    "00:00",   "07:30", "09:00",   "15:00", "16:30",   "24:00", FALSE,
  0.8,  "09:30",   "15:00", "16:30",   "20:30", NA,        NA,      FALSE,
  0.7,  "16:30",   "19:30", NA,        NA,      NA,        NA,      FALSE,
  0.6,  "16:30",   "24:00", NA,        NA,      NA,        NA,      TRUE,    # continues 00:00–07:30 next day
  0.5,  "09:30",   "12:00", "16:30",   "19:00", NA,        NA,      FALSE,
  0.4,  "09:30",   "15:00", NA,        NA,      NA,        NA,      FALSE,
  0.3,  "09:30",   "12:30", NA,        NA,      NA,        NA,      FALSE,
  0.2,  "09:00",   "11:00", NA,        NA,      NA,        NA,      FALSE
)


# =============================================================================
# SECTION 1 — Site event date lists
# =============================================================================
# The `management` list stores the dates on which each event type occurred
# at each site and year.  Only the dates are stored here; time windows are
# resolved in Sections 2 and 3 using the code_map above.
#
# Structure:
#   management[["<SITE>"]]$events[["<YEAR>"]][["<EventType>"]] <- c("YYYY-MM-DD", ...)
#
# Supported event types: "Grazing", "Fertiliser", "Slurry", "Harvest"
# Add or remove event types as needed for your site.

management <- list(

  # ---------------------------------------------------------------------------
  # SITE 1  (replace "SITE1" with your site identifier)
  # ---------------------------------------------------------------------------
  "SITE1" = list(
    tz = "UTC",
    events = list(
      "2023" = list(
        Grazing     = c(
          # Add one "YYYY-MM-DD" string per grazing day
          "YYYY-MM-DD"   # <-- replace with your first grazing date
          # "YYYY-MM-DD" # <-- add more rows as needed
        ),
        Fertiliser  = c(
          # Add one "YYYY-MM-DD" string per fertiliser application date
          "YYYY-MM-DD"   # <-- replace
        )
      ),
      "2024" = list(
        Grazing     = c(
          "YYYY-MM-DD"   # <-- replace
        ),
        Fertiliser  = c(
          "YYYY-MM-DD"   # <-- replace
        )
      )
    )
  ),

  # ---------------------------------------------------------------------------
  # SITE 2  (replace "SITE2" with your site identifier, copy block for more)
  # ---------------------------------------------------------------------------
  "SITE2" = list(
    tz = "UTC",
    events = list(
      "2023" = list(
        Grazing    = c(
          "YYYY-MM-DD"   # <-- replace
        ),
        Fertiliser = c(
          "YYYY-MM-DD"   # <-- replace
        )
      ),
      "2024" = list(
        Grazing    = c(
          "YYYY-MM-DD"   # <-- replace
        ),
        Fertiliser = c(
          "YYYY-MM-DD"   # <-- replace
        ),
        Slurry     = c(
          # Slurry is treated separately from mineral fertiliser (excluded from N step function)
          "YYYY-MM-DD"   # <-- replace, or remove this block if no slurry events
        )
      )
    )
  )
)

# De-duplicate grazing dates within each year (in case the same date appears
# more than once due to data entry errors or overlapping records).
management[["SITE1"]]$events[["2023"]][["Grazing"]] <- unique(as.Date(management[["SITE1"]]$events[["2023"]][["Grazing"]]))
management[["SITE1"]]$events[["2024"]][["Grazing"]] <- unique(as.Date(management[["SITE1"]]$events[["2024"]][["Grazing"]]))
management[["SITE2"]]$events[["2023"]][["Grazing"]] <- unique(as.Date(management[["SITE2"]]$events[["2023"]][["Grazing"]]))
management[["SITE2"]]$events[["2024"]][["Grazing"]] <- unique(as.Date(management[["SITE2"]]$events[["2024"]][["Grazing"]]))


# =============================================================================
# SECTION 2 — Per-site, per-year grazing detail tables
# =============================================================================
# For each grazing day, record which intensity code was observed.  The code
# is looked up in code_map (Section 0) to produce the exact half-hour windows.
#
# If you have an explicit time window (e.g. "09:00-12:00" from field notes)
# rather than a code, you can use that directly — rows with a t_start/t_end
# pair bypass the code_map lookup.
#
# Column descriptions:
#   date_chr  — date as a string in DD.MM.YYYY format
#   detail    — either "code X.X" (e.g. "code 1") or "HH:MM-HH:MM"
#
# After the tribble, date and code are parsed automatically.

# --- SITE1 — Year 2023 -------------------------------------------------------
site1_grazing_detail_2023 <- tibble::tribble(
  ~date_chr,      ~detail,
  # Each row = one grazing day.  "code X.X" maps to the code_map above.
  # "DD.MM.YYYY",  "code 1",     # full-day grazing (highest intensity)
  # "DD.MM.YYYY",  "code 0.8",   # morning + afternoon grazing
  # "DD.MM.YYYY",  "code 0.5",   # mid-morning + early evening
  # "DD.MM.YYYY",  "code 0.3",   # short morning window only
  # "DD.MM.YYYY",  "09:30-11:00",# explicit window (overrides code_map)
  "DD.MM.YYYY",  "code 0.0"    # <-- replace with your first grazing record
) %>%
  mutate(date = lubridate::dmy(date_chr),
         code = readr::parse_number(detail)) %>%
  dplyr::select(date, code)

# --- SITE1 — Year 2024 -------------------------------------------------------
# For years with a mix of intensity codes AND explicit time windows, add
# t_start and t_end columns as shown in this template.
site1_grazing_detail_2024 <- tibble::tribble(
  ~date_chr,      ~detail,
  # "DD.MM.YYYY",  "code 1",
  # "DD.MM.YYYY",  "09:00-12:00",  # explicit morning window
  "DD.MM.YYYY",  "code 0.0"    # <-- replace
) %>%
  mutate(
    date    = lubridate::dmy(date_chr),
    code    = readr::parse_number(detail),
    # Extract explicit HH:MM-HH:MM time strings if present
    times   = stringr::str_extract(detail, "\\b\\d{1,2}:\\d{2}-\\d{1,2}:\\d{2}\\b"),
    times   = stringr::str_replace_all(times, "^(\\d):(\\d{2})", "0\\1:\\2"),
    times   = stringr::str_replace_all(times, "-(\\d):(\\d{2})$", "-0\\1:\\2"),
    t_start = dplyr::if_else(!is.na(times), stringr::str_sub(times, 1, 5), NA_character_),
    t_end   = dplyr::if_else(!is.na(times), stringr::str_sub(times, 7, 11), NA_character_)
  ) %>%
  dplyr::select(date, code, t_start, t_end)

# --- SITE2 — Year 2023 -------------------------------------------------------
site2_grazing_detail_2023 <- tibble::tribble(
  ~date_chr,      ~detail,
  "DD.MM.YYYY",  "code 0.0"    # <-- replace
) %>%
  mutate(date = lubridate::dmy(date_chr),
         code = readr::parse_number(detail)) %>%
  dplyr::select(date, code)

# --- SITE2 — Year 2024 -------------------------------------------------------
site2_grazing_detail_2024 <- tibble::tribble(
  ~date_chr,      ~detail,
  "DD.MM.YYYY",  "code 0.0"    # <-- replace
) %>%
  mutate(
    date    = lubridate::dmy(date_chr),
    code    = readr::parse_number(detail),
    times   = stringr::str_extract(detail, "\\b\\d{2}:\\d{2}-\\d{2}:\\d{2}\\b"),
    t_start = dplyr::if_else(!is.na(times), stringr::str_sub(times, 1, 5), NA_character_),
    t_end   = dplyr::if_else(!is.na(times), stringr::str_sub(times, 7, 11), NA_character_)
  ) %>%
  dplyr::select(date, code, t_start, t_end)


# =============================================================================
# SECTION 3 — Helper functions  (DO NOT EDIT)
# =============================================================================
# These functions convert the date-and-code tables above into standardised
# event tables with exact start_time and end_time columns.
# They are data-independent — edit only Section 0, 1, and 2 above.

# Pad HH:MM strings to always have a leading zero (e.g. "9:30" -> "09:30")
pad_times <- function(df) {
  df %>%
    mutate(
      start_time = stringr::str_pad(start_time, 5, pad = "0"),
      end_time   = stringr::str_pad(end_time,   5, pad = "0")
    )
}

# Expand per-day intensity codes into one-to-three time windows per day.
# carry_next = TRUE (code 0.6) adds a 00:00–07:30 window on the following day.
# allowed_next_day restricts carry-over to dates that are known grazing days.
expand_code_windows_grouped <- function(df_code, allowed_next_day = NULL,
                                         event_label = "Grazing") {
  dfc <- df_code %>%
    dplyr::left_join(code_map, by = "code") %>%
    dplyr::mutate(origin_date = date,
                  occ_group   = paste0("occ_", format(origin_date, "%Y-%m-%d")))

  w <- dplyr::bind_rows(
    dfc %>% dplyr::filter(!is.na(t1_start)) %>%
      dplyr::transmute(date, event = event_label,
                       start_time = t1_start, end_time = t1_end, occ_group),
    dfc %>% dplyr::filter(!is.na(t2_start)) %>%
      dplyr::transmute(date, event = event_label,
                       start_time = t2_start, end_time = t2_end, occ_group),
    dfc %>% dplyr::filter(!is.na(t3_start)) %>%
      dplyr::transmute(date, event = event_label,
                       start_time = t3_start, end_time = t3_end, occ_group)
  )

  # Carry-over window for code 0.6 (same occurrence group as the preceding day)
  carry <- dfc %>%
    dplyr::filter(carry_next) %>%
    dplyr::transmute(
      date       = origin_date + lubridate::days(1),
      event      = event_label,
      start_time = "00:00",
      end_time   = "07:30",
      occ_group
    )

  if (!is.null(allowed_next_day))
    carry <- carry %>% dplyr::filter(date %in% allowed_next_day)

  dplyr::bind_rows(w, carry) %>% pad_times()
}

# Number grazing occurrences chronologically within each event category.
# Occurrences are identified by occ_group (same-day windows share one group).
number_events_all <- function(ev_tbl, tz = "UTC") {
  if (!nrow(ev_tbl)) return(ev_tbl)

  ev_tbl <- ev_tbl %>%
    pad_times() %>%
    mutate(
      event    = as.character(event),
      category = stringr::str_remove(event, "_\\d+$"),
      has_idx  = stringr::str_detect(event, "_\\d+$")
    )

  if (!"occ_group" %in% names(ev_tbl)) ev_tbl$occ_group <- NA_character_

  ev_tbl <- ev_tbl %>%
    mutate(
      occ_group = dplyr::case_when(
        has_idx   ~ paste0(category, "::", event),
        is.na(occ_group) | occ_group == "" ~
          paste0(category, "::", format(as.Date(date), "%Y-%m-%d")),
        TRUE ~ occ_group
      ),
      start_ct = as.POSIXct(paste(date, start_time), tz = tz)
    )

  grp_order <- ev_tbl %>%
    group_by(category, occ_group) %>%
    summarise(group_start = min(start_ct, na.rm = TRUE), .groups = "drop") %>%
    arrange(category, group_start) %>%
    group_by(category) %>%
    mutate(idx = dplyr::row_number()) %>%
    ungroup()

  ev_tbl %>%
    left_join(grp_order, by = c("category", "occ_group")) %>%
    mutate(event = if_else(has_idx, event, paste0(category, "_", idx))) %>%
    dplyr::select(date, event, start_time, end_time) %>%
    arrange(as.Date(date), event, start_time)
}

# Sort event table and remove exact duplicates (same date, event, start, end).
finalize_site_events <- function(ev_tbl) {
  ev_tbl %>%
    dplyr::mutate(
      date       = as.Date(date),
      start_time = stringr::str_pad(start_time, 5, pad = "0"),
      end_time   = stringr::str_pad(end_time,   5, pad = "0")
    ) %>%
    dplyr::arrange(date, event, start_time, end_time) %>%
    dplyr::distinct(date, event, start_time, end_time, .keep_all = TRUE)
}


# =============================================================================
# SECTION 4 — Build site-year event tables  (DO NOT EDIT LOGIC)
# =============================================================================
# This section calls the helper functions above to expand the date lists and
# detail tables into the final event table format expected by the pipeline.
# If you add a new site, copy one of the blocks below, update the variable
# names, and add a finalize_site_events() call at the end.

# --- SITE1 — 2023 ------------------------------------------------------------
allowed_site1_2023 <- as.Date(management[["SITE1"]]$events[["2023"]][["Grazing"]])

site1_grazing_windows_2023 <- expand_code_windows_grouped(
  site1_grazing_detail_2023,
  allowed_next_day = allowed_site1_2023
)

site1_fertiliser_windows_2023 <- tibble::tibble(
  date       = as.Date(management[["SITE1"]]$events[["2023"]][["Fertiliser"]]),
  event      = "Fertiliser",
  start_time = "10:30",
  end_time   = "11:00"
) %>% pad_times()

site1_events_2023 <- dplyr::bind_rows(
  site1_grazing_windows_2023,
  site1_fertiliser_windows_2023
) %>% number_events_all()

# --- SITE1 — 2024 ------------------------------------------------------------
# Rows with an explicit t_start keep their own occurrence group.
gw_explicit_site1_2024 <- site1_grazing_detail_2024 %>%
  dplyr::filter(!is.na(t_start)) %>%
  dplyr::transmute(
    date, event = "Grazing",
    start_time = t_start, end_time = t_end,
    occ_group  = paste0("occ_", format(date, "%Y-%m-%d"))
  ) %>% pad_times()

gw_code_base_site1_2024 <- site1_grazing_detail_2024 %>%
  dplyr::filter(is.na(t_start) & !is.na(code)) %>%
  dplyr::select(date, code)

allowed_site1_2024 <- as.Date(management[["SITE1"]]$events[["2024"]][["Grazing"]])
gw_code_site1_2024 <- expand_code_windows_grouped(
  gw_code_base_site1_2024,
  allowed_next_day = allowed_site1_2024
)

site1_fertiliser_windows_2024 <- tibble::tibble(
  date       = as.Date(management[["SITE1"]]$events[["2024"]][["Fertiliser"]]),
  event      = "Fertiliser",
  start_time = "10:30",
  end_time   = "11:00"
) %>% pad_times()

site1_events_2024 <- dplyr::bind_rows(
  gw_explicit_site1_2024,
  gw_code_site1_2024,
  site1_fertiliser_windows_2024
) %>% number_events_all()

# --- SITE2 — 2023 ------------------------------------------------------------
allowed_site2_2023 <- as.Date(management[["SITE2"]]$events[["2023"]][["Grazing"]])

site2_grazing_windows_2023 <- expand_code_windows_grouped(
  site2_grazing_detail_2023,
  allowed_next_day = allowed_site2_2023
)

site2_fertiliser_windows_2023 <- tibble::tibble(
  date       = as.Date(management[["SITE2"]]$events[["2023"]][["Fertiliser"]]),
  event      = "Fertiliser",
  start_time = "10:30",
  end_time   = "11:00"
) %>% pad_times()

site2_events_2023 <- dplyr::bind_rows(
  site2_grazing_windows_2023,
  site2_fertiliser_windows_2023
) %>% number_events_all()

# --- SITE2 — 2024 ------------------------------------------------------------
gw_explicit_site2_2024 <- site2_grazing_detail_2024 %>%
  dplyr::filter(!is.na(t_start)) %>%
  dplyr::transmute(
    date, event = "Grazing",
    start_time = t_start, end_time = t_end,
    occ_group  = paste0("occ_", format(date, "%Y-%m-%d"))
  ) %>% pad_times()

gw_code_base_site2_2024 <- site2_grazing_detail_2024 %>%
  dplyr::filter(is.na(t_start) & !is.na(code)) %>%
  dplyr::select(date, code)

allowed_site2_2024 <- as.Date(management[["SITE2"]]$events[["2024"]][["Grazing"]])
gw_code_site2_2024 <- expand_code_windows_grouped(
  gw_code_base_site2_2024,
  allowed_next_day = allowed_site2_2024
)

site2_grazing_windows_2024 <- dplyr::bind_rows(gw_explicit_site2_2024, gw_code_site2_2024)

# Separate mineral fertiliser from slurry (slurry excluded from N step function)
site2_fert_all_2024   <- as.Date(management[["SITE2"]]$events[["2024"]][["Fertiliser"]])
site2_slurry_2024     <- as.Date(management[["SITE2"]]$events[["2024"]][["Slurry"]])
site2_mineral_fert_24 <- setdiff(site2_fert_all_2024, site2_slurry_2024)

site2_fertiliser_windows_2024 <- tibble::tibble(
  date       = site2_mineral_fert_24,
  event      = "Fertiliser",
  start_time = "10:30",
  end_time   = "11:00"
) %>% pad_times()

site2_slurry_windows_2024 <- tibble::tibble(
  date       = site2_slurry_2024,
  event      = "Slurry",
  start_time = "10:30",
  end_time   = "11:00"
) %>% pad_times()

site2_events_2024 <- dplyr::bind_rows(
  site2_grazing_windows_2024,
  site2_fertiliser_windows_2024,
  site2_slurry_windows_2024
) %>% number_events_all()


# --- Finalise all site-year event tables -------------------------------------
site1_events_2023 <- finalize_site_events(site1_events_2023)
site1_events_2024 <- finalize_site_events(site1_events_2024)
site2_events_2023 <- finalize_site_events(site2_events_2023)
site2_events_2024 <- finalize_site_events(site2_events_2024)

# ======================= end ManagementEvents.R ==============================
