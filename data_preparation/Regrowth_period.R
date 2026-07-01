# Regrowth_period.R — regrowth-period, grazing-group, cow and N annotations
#
# Site-level utility (run once per site). Derives grazing "groups" and the
# regrowth window after each group from the management event tables, then adds
# per-group dummies, the cow count during grazing, and the N step function to a
# gap-filling data frame, and saves the annotated frame back to disk.
#
# NOTE: the source() and readRDS()/saveRDS() paths below use an older repository
# layout ("03. Data_Preparation/...", "01. Data_Main/...") and will not resolve
# against the current structure — update them to your paths before running.
#
# To switch site, edit the SITE-LEVEL block (site_id, rds_name, events_*,
# cow_windows, N_windows) from the jc1_* selections to the jc2_* ones.

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

# hard-coded management event tables (jc1_events_*, jc2_events_*)
source(here::here("03. Data_Preparation/01d. ManagementEvents.R"))

# site-level cow-count windows and N schedules (JC1 and JC2 both provided)
cow_windows_jc1 <- tibble::tribble(
  ~start_date,   ~end_date,     ~cow_number,
  # ---- 2020 ----
  "2020-02-04",  "2020-02-07",  30L,   # ramp 27 -> 32; start value
  "2020-02-10",  "2020-02-10",  37L,
  "2020-03-03",  "2020-03-07",  21L,
  "2020-03-09",  "2020-03-09",  21L,
  "2020-03-12",  "2020-03-13",  21L,
  "2020-03-16",  "2020-03-17",  21L,
  "2020-03-19",  "2020-03-22",  22L,   # tied with 24; start value
  "2020-04-10",  "2020-04-18",  29L,
  "2020-05-03",  "2020-05-09",  29L,
  "2020-05-25",  "2020-06-02",  29L,
  "2020-06-17",  "2020-06-25",  27L,
  "2020-07-09",  "2020-07-18",  27L,
  "2020-08-01",  "2020-08-12",  22L,
  "2020-09-01",  "2020-09-21",  19L,   # tied with 20; start value
  "2020-10-21",  "2020-10-30",  29L,
  "2020-11-02",  "2020-11-07",  31L,
  # ---- 2023 ----
  "2023-02-19",  "2023-02-28",  12L,
  "2023-03-01",  "2023-03-11",  11L,
  "2023-04-10",  "2023-04-27",  17L,
  "2023-05-11",  "2023-05-13",  19L,
  "2023-05-16",  "2023-05-21",  19L,
  "2023-05-29",  "2023-05-31",  19L,
  "2023-06-09",  "2023-06-12",  19L,
  "2023-06-19",  "2023-06-30",  19L,
  "2023-07-09",  "2023-07-14",  19L,
  "2023-08-16",  "2023-08-19",  29L,
  "2023-08-22",  "2023-09-08",  14L,
  "2023-10-16",  "2023-10-16",  77L,
  # ---- 2024 ----
  "2024-02-27",  "2024-02-27",  18L,
  "2024-02-29",  "2024-02-29",  18L,
  "2024-03-03",  "2024-03-12",  18L,
  "2024-04-16",  "2024-04-30",  18L,
  "2024-05-01",  "2024-05-01",  14L,
  "2024-05-06",  "2024-05-08",  14L,
  "2024-05-10",  "2024-05-10",  14L,
  "2024-05-17",  "2024-05-28",  21L,
  "2024-06-09",  "2024-06-27",  21L,
  "2024-07-23",  "2024-08-09",  11L
)

N_windows_jc1 <- tibble::tribble(
  ~fert_date,    ~N_kg_ha,
  # ---- 2020 (CAN only) ----
  "02/04/2020",   50,
  "11/05/2020",   40,
  "03/06/2020",   27,
  "29/06/2020",   20,
  "14/08/2020",   27,
  "14/09/2020",   27,
  # ---- 2023 ----
  "20/02/2023",   18.5,
  # "19/05/2023" Slurry excluded on purpose
  # ---- 2024 ----
  "02/04/2024",   30.4,
  "07/05/2024",   23,
  "29/05/2024",   17,
  "20/06/2024",   10,
  "06/08/2024",   10
  # "22/08/2024" Slurry excluded on purpose
)

cow_windows_jc2 <- tibble::tribble(
  ~start_date, ~end_date, ~cow_number,
  "2023-03-21", "2023-03-24", 43L,
  "2023-04-24", "2023-04-25", 43L,
  "2023-05-21", "2023-05-24", 42L,
  "2023-06-21", "2023-06-25", 42L,
  "2023-07-21", "2023-07-22", 84L,
  "2023-08-17", "2023-08-21", 74L,
  "2023-09-23", "2023-09-23", 97L,
  "2023-09-25", "2023-09-26", 97L,
  "2023-09-28", "2023-09-29", 97L,
  "2023-10-03", "2023-10-04", 97L,
  "2024-02-07", "2024-02-07", 88L,
  "2024-02-12", "2024-02-13", 88L,
  "2024-04-15", "2024-04-18", 44L,
  "2024-05-22", "2024-05-28", 42L,
  "2024-06-20", "2024-06-25", 42L,
  "2024-07-26", "2024-07-30", 41L
)

N_windows_jc2 <- tibble::tribble(
  ~fert_date,   ~N_kg_ha,
  "20/02/2023",  30,
  "04/04/2023",  35,
  "26/04/2023",  32,
  "22/05/2023",  30,
  "20/06/2023",  30,
  "12/07/2023",  27,
  "24/08/2023",  27,
  "08/09/2023",  21,
  "02/04/2024",  38,
  # "19/04/2024" Slurry excluded on purpose
  "07/05/2024",  30.4,
  "29/05/2024",  25,
  "19/06/2024",  28,
  "10/07/2024",  20,
  "06/08/2024",  25
)
# choose the site here (swap the jc1_* selections for jc2_* to run JC2)
site_id     <- "JC1"
rds_name    <- "JC1_2023_2024_df.rds"
events_2023 <- jc1_events_2023
events_2024 <- jc1_events_2024
cow_windows <- cow_windows_jc1
N_windows   <- N_windows_jc1

# 1) grazing groups: bind both years and cluster events into groups
#    separated by gaps of more than 10 days
grazing_events <- bind_rows(
  events_2023,
  events_2024
) %>%
  filter(str_detect(event, "^Grazing")) %>%
  mutate(
    date_dt = as.POSIXct(date, tz = "UTC"),
    start_dt = date_dt +
      hm(if_else(start_time == "24:00", "00:00", start_time)) +
      days(if_else(start_time == "24:00", 1, 0)),
    end_dt = date_dt +
      hm(if_else(end_time == "24:00", "00:00", end_time)) +
      days(if_else(end_time == "24:00", 1, 0))
  ) %>%
  arrange(start_dt) %>%
  mutate(
    grazing_group = cumsum(
      if_else(
        is.na(lag(start_dt)) |
          as.numeric(difftime(start_dt, lag(start_dt), units = "days")) > 10,
        1L, 0L
      )
    )
  )

grazing_groups <- grazing_events %>%
  group_by(grazing_group) %>%
  summarise(
    group_start    = min(start_dt),
    last_event_end = max(end_dt),
    .groups = "drop"
  ) %>%
  arrange(group_start) %>%
  mutate(
    group_end = coalesce(lead(group_start), last_event_end)
  )

# 2) read the gap-filling data frame and add group dummies, cows and N
df <- readRDS(here::here("01. Data_Main/Data_GapFilling", site_id, rds_name))

tz_df <- attr(df$timestamp, "tzone")
if (is.null(tz_df) || tz_df == "") tz_df <- "UTC"

grazing_groups_tz <- grazing_groups %>%
  mutate(
    grazing_group  = as.integer(grazing_group),
    group_start    = with_tz(group_start,    tz_df),
    last_event_end = with_tz(last_event_end, tz_df)
  ) %>%
  arrange(group_start)

first_start <- grazing_groups_tz$group_start[1]
last_end    <- grazing_groups_tz$last_event_end[nrow(grazing_groups_tz)]

df <- df %>%
  mutate(
    pre_grazing  = as.integer(timestamp < first_start),
    post_grazing = as.integer(timestamp >= last_end),
    .row_id = row_number()
  )

# grazing_group_i dummy: 1 within [group_start, last_event_end)
grazing_group_dummies <- df %>%
  select(.row_id, timestamp) %>%
  crossing(grazing_groups_tz %>% select(grazing_group, group_start, last_event_end)) %>%
  mutate(
    val = as.integer(timestamp >= group_start & timestamp < last_event_end),
    var = paste0("grazing_group_", grazing_group)
  ) %>%
  select(.row_id, var, val) %>%
  pivot_wider(names_from = var, values_from = val, values_fill = 0)

# growth_i dummy: 1 within [last_event_end_i, next group_start)
growth_def <- grazing_groups_tz %>%
  transmute(
    grazing_group,
    growth_start = last_event_end,
    growth_end   = lead(group_start)
  ) %>%
  mutate(growth_end = coalesce(growth_end, growth_start))

growth_dummies <- df %>%
  select(.row_id, timestamp) %>%
  crossing(growth_def) %>%
  mutate(
    val = as.integer(timestamp >= growth_start & timestamp < growth_end),
    var = paste0("growth_", grazing_group)
  ) %>%
  select(.row_id, var, val) %>%
  pivot_wider(names_from = var, values_from = val, values_fill = 0)

df <- df %>%
  left_join(grazing_group_dummies, by = ".row_id") %>%
  left_join(growth_dummies,        by = ".row_id")

# cow count: value active only where Grazing == 1
cow_windows2 <- cow_windows %>%
  mutate(
    start_date = ymd(start_date),
    end_date   = ymd(end_date),
    start_dt   = as.POSIXct(start_date, tz = tz_df),
    end_dt     = as.POSIXct(end_date + days(1), tz = tz_df)
  )

df <- df %>%
  left_join(
    df %>%
      select(.row_id, timestamp, Grazing) %>%
      crossing(cow_windows2 %>% select(start_dt, end_dt, cow_number)) %>%
      mutate(hit = Grazing == 1 & timestamp >= start_dt & timestamp < end_dt) %>%
      group_by(.row_id) %>%
      summarise(
        cows = if_else(any(hit), first(cow_number[hit]), 0L),
        .groups = "drop"
      ),
    by = ".row_id"
  )

# N: step function that changes at each fertilisation event

N_windows2 <- N_windows %>%
  mutate(
    fert_date = str_replace_all(fert_date, "\\.", "/"),
    fert_date = dmy(fert_date)
  )

# use the *observed* first half-hour where Fertiliser==1 on each fertilisation date (fallback to midnight)
fert_starts <- df %>%
  filter(Fertiliser == 1) %>%
  mutate(fert_date = as.Date(timestamp)) %>%
  group_by(fert_date) %>%
  summarise(start_dt_obs = min(timestamp), .groups = "drop")

fert_schedule <- N_windows2 %>%
  left_join(fert_starts, by = "fert_date") %>%
  mutate(
    start_dt = coalesce(start_dt_obs, as.POSIXct(fert_date, tz = tz_df))
  ) %>%
  arrange(start_dt)

n_starts <- fert_schedule$start_dt
n_vals   <- fert_schedule$N_kg_ha

df <- df %>%
  mutate(
    .n_idx  = findInterval(timestamp, n_starts),      # 0 before first start
    .n_idx2 = if_else(.n_idx == 0L, 1L, .n_idx),      # safe for indexing
    N = if_else(.n_idx == 0L, 0, as.numeric(n_vals[.n_idx2]))
  ) %>%
  select(-.n_idx, -.n_idx2)


# save the annotated data frame back to the same path
out_path <- here::here("01. Data_Main/Data_GapFilling", site_id, rds_name)
saveRDS(df, out_path)
