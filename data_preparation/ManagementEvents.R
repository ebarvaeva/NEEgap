# ManagementEvents.R — hard-coded JC1 & JC2 grazing / fertiliser event tables
#
# Builds the per-site, per-year event tables (jc1_events_2020/2023/2024,
# jc2_events_2023/2024) that generate_management_csvs.R exports to CSV. Grazing
# "codes" map to one or more half-hour windows per day via code_map; fertiliser
# events use a fixed 10:30-11:00 window. This script defines tables in memory
# and saves nothing itself.

library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(readr)
library(lubridate)
library(purrr)

# grazing "code" -> up to three daily windows (24:00 is normalised below)
code_map <- tibble::tribble(
  ~code, ~t1_start, ~t1_end, ~t2_start, ~t2_end, ~t3_start, ~t3_end, ~carry_next,
  1,    "00:00",   "07:30", "09:00",   "15:00", "16:30",   "24:00", FALSE,
  0.8,  "09:30",   "15:00", "16:30",   "20:30", NA,        NA,      FALSE,
  0.7,  "16:30",   "19:30", NA,        NA,      NA,        NA,      FALSE,
  0.6,  "16:30",   "24:00", NA,        NA,      NA,        NA,      TRUE,   # plus next-day 00:00–07:30
  0.5,  "09:30",   "12:00", "16:30",   "19:00", NA,        NA,      FALSE,
  0.4,  "09:30",   "15:00", NA,        NA,      NA,        NA,      FALSE,
  0.3,  "09:30",   "12:30", NA,        NA,      NA,        NA,      FALSE,
  0.2,  "09:00",   "11:00", NA,        NA,      NA,        NA,      FALSE
)

# site-specific grazing / fertiliser event dates
management <- list(
  "JC1" = list(tz = "UTC",
               events = list(
                 "2020" = list(
                   Grazing = c(
                     "2020-02-04","2020-02-05","2020-02-06","2020-02-07","2020-02-10",
                     "2020-03-03","2020-03-04","2020-03-05","2020-03-06","2020-03-07",
                     "2020-03-09","2020-03-12","2020-03-13","2020-03-16","2020-03-17",
                     "2020-03-19","2020-03-20","2020-03-21","2020-03-22","2020-04-10",
                     "2020-04-11","2020-04-12","2020-04-13","2020-04-14","2020-04-15",
                     "2020-04-16","2020-04-17","2020-04-18","2020-05-03","2020-05-04",
                     "2020-05-05","2020-05-06","2020-05-07","2020-05-08","2020-05-09",
                     "2020-05-25","2020-05-26","2020-05-27","2020-05-28","2020-05-29",
                     "2020-05-30","2020-05-31","2020-06-01","2020-06-02","2020-06-17",
                     "2020-06-18","2020-06-19","2020-06-20","2020-06-21","2020-06-22",
                     "2020-06-23","2020-06-24","2020-06-25","2020-07-09","2020-07-10",
                     "2020-07-11","2020-07-12","2020-07-13","2020-07-14","2020-07-15",
                     "2020-07-16","2020-07-17","2020-07-18","2020-08-01","2020-08-02",
                     "2020-08-03","2020-08-04","2020-08-05","2020-08-06","2020-08-07",
                     "2020-08-08","2020-08-09","2020-08-10","2020-08-11","2020-08-12",
                     "2020-09-01","2020-09-02","2020-09-03","2020-09-04","2020-09-05",
                     "2020-09-06","2020-09-07","2020-09-08","2020-09-09","2020-09-10",
                     "2020-09-11","2020-09-12","2020-09-13","2020-09-14","2020-09-15",
                     "2020-09-16","2020-09-17","2020-09-18","2020-09-19","2020-09-20",
                     "2020-09-21","2020-10-21","2020-10-22","2020-10-23","2020-10-24",
                     "2020-10-25","2020-10-26","2020-10-27","2020-10-28","2020-10-29",
                     "2020-10-30","2020-11-02","2020-11-03","2020-11-04","2020-11-05",
                     "2020-11-06","2020-11-07"
                   ),
                   Fertiliser = c("2020-03-03","2020-04-02","2020-05-04","2020-05-11","2020-05-25",
                                  "2020-06-03","2020-06-29","2020-08-14","2020-09-01","2020-09-14")
                 ),
                 "2023" = list(
                   Grazing = c(
                     "2023-02-19","2023-02-20","2023-02-21","2023-02-22","2023-02-23","2023-02-24",
                     "2023-02-25","2023-02-26","2023-02-27","2023-02-28","2023-03-01","2023-03-02",
                     "2023-03-03","2023-03-04","2023-03-05","2023-03-06","2023-03-07","2023-03-08",
                     "2023-03-09","2023-03-10","2023-03-11","2023-04-10","2023-04-11",
                     "2023-04-14","2023-04-15","2023-04-16","2023-04-17","2023-04-18",
                     "2023-04-19","2023-04-20","2023-04-21","2023-04-22","2023-04-23","2023-04-24",
                     "2023-04-25","2023-04-26","2023-04-27","2023-05-11","2023-05-12","2023-05-13",
                     "2023-05-16","2023-05-17","2023-05-18","2023-05-19","2023-05-20","2023-05-21",
                     "2023-05-22",
                     "2023-05-29","2023-05-30","2023-05-31","2023-06-09","2023-06-10","2023-06-11",
                     "2023-06-12","2023-06-19","2023-06-20","2023-06-21","2023-06-22","2023-06-23",
                     "2023-06-24","2023-06-25","2023-06-26","2023-06-27","2023-06-28","2023-06-29",
                     "2023-06-30","2023-07-09","2023-07-10","2023-07-11","2023-07-12","2023-07-13",
                     "2023-07-14","2023-07-20","2023-07-21","2023-07-22","2023-07-23","2023-07-24",
                     "2023-07-25","2023-07-26","2023-07-27","2023-07-28","2023-07-29","2023-07-30",
                     "2023-07-31","2023-08-01","2023-08-16","2023-08-17","2023-08-18","2023-08-19",
                     "2023-08-22","2023-08-23","2023-08-24","2023-08-25","2023-08-26","2023-08-27",
                     "2023-08-28","2023-08-29","2023-08-30","2023-08-31","2023-09-01","2023-09-02",
                     "2023-09-03","2023-09-04","2023-09-05","2023-09-06","2023-09-07","2023-09-08",
                     "2023-10-16"
                   ),
                   Fertiliser = c("2023-02-20","2023-02-28","2023-05-03","2023-05-24","2023-05-19")
                 ),
                 "2024" = list(
                   Grazing = c(
                     "2024-02-27","2024-02-29","2024-02-29","2024-03-03","2024-03-04","2024-03-05",
                     "2024-03-06","2024-03-07","2024-03-08","2024-03-09","2024-03-11","2024-03-12",
                     "2024-04-16","2024-04-17","2024-04-18","2024-04-19","2024-04-20","2024-04-21",
                     "2024-04-22","2024-04-23","2024-04-24","2024-04-25","2024-04-26","2024-04-27",
                     "2024-04-28","2024-04-29","2024-04-30","2024-05-01","2024-05-06","2024-05-07",
                     "2024-05-08","2024-05-10","2024-05-17","2024-05-18","2024-05-19","2024-05-20",
                     "2024-05-21","2024-05-22","2024-05-23","2024-05-24","2024-05-25","2024-05-26",
                     "2024-05-27","2024-05-28","2024-06-09","2024-06-10","2024-06-11","2024-06-12",
                     "2024-06-13","2024-06-14","2024-06-15","2024-06-16","2024-06-17","2024-06-18",
                     "2024-06-19","2024-06-20","2024-06-21","2024-06-22","2024-06-23","2024-06-24",
                     "2024-06-25","2024-06-26","2024-06-27","2024-07-23","2024-07-24","2024-07-25",
                     "2024-07-26","2024-07-27","2024-07-28","2024-07-29","2024-07-30","2024-07-31",
                     "2024-08-01","2024-08-02","2024-08-03","2024-08-04","2024-08-05","2024-08-06",
                     "2024-08-07","2024-08-08","2024-08-09","2024-09-25","2024-09-26","2024-09-27",
                     "2024-09-28","2024-10-12","2024-10-13"
                   ),
                   Fertiliser = c("2024-04-02","2024-05-07","2024-05-29","2024-06-20","2024-08-06","2024-09-10")
                 )
               )),
  "JC2" = list(tz = "UTC",
               events = list(
                 "2023" = list(
                   Grazing = c(
                     "2023-03-21","2023-03-22","2023-03-23","2023-03-24","2023-04-24","2023-04-25",
                     "2023-05-21","2023-05-22","2023-05-23","2023-05-24","2023-06-21",
                     "2023-06-22","2023-06-23","2023-06-24","2023-06-25","2023-07-21",
                     "2023-07-22","2023-08-17","2023-08-18","2023-08-19","2023-08-20",
                     "2023-08-21","2023-09-23","2023-09-24","2023-09-25","2023-09-26",
                     "2023-09-27","2023-09-28",
                     "2023-09-29","2023-10-03","2023-10-04"
                   ),
                   Fertiliser = c(
                     "2023-02-20","2023-02-28","2023-04-04","2023-04-26","2023-05-03","2023-05-22","2023-06-20",
                     "2023-07-12","2023-08-24","2023-09-08", "2023-11-30"
                   )
                 ),
                 "2024" = list(
                   Grazing = c(
                     "2024-02-07","2024-02-12","2024-02-13","2024-04-15","2024-04-16","2024-04-17",
                     "2024-04-18","2024-05-22","2024-05-23","2024-05-24","2024-05-25","2024-05-27",
                     "2024-05-27","2024-05-28","2024-06-20","2024-06-21","2024-06-22","2024-06-23",
                     "2024-06-24","2024-06-25","2024-07-26","2024-07-27","2024-07-28","2024-07-29",
                     "2024-07-30","2024-09-09","2024-09-10","2024-09-11","2024-09-12","2024-09-13",
                     "2024-09-14","2024-09-15"
                   ),
                   Fertiliser = c(
                     "2024-04-02","2024-04-19","2024-05-07","2024-05-29","2024-06-19",
                     "2024-07-10","2024-08-06","2024-09-10"
                   ),
                   Slurry = c("2024-04-19")
                 )
               ))
)

# JC1 2023 grazing detail: date -> grazing code
jc1_grazing_detail_2023 <- tibble::tribble(
  ~date_chr,    ~detail,
  "19.02.2023", "code 1",
  "20.02.2023", "code 1",
  "21.02.2023", "code 1",
  "22.02.2023", "code 1",
  "23.02.2023", "code 0.6",
  "24.02.2023", "code 0.6",
  "25.02.2023", "code 1",
  "26.02.2023", "code 1",
  "27.02.2023", "code 1",
  "28.02.2023", "code 1",
  "01.03.2023", "code 1",
  "02.03.2023", "code 1",
  "03.03.2023", "code 1",
  "04.03.2023", "code 1",
  "05.03.2023", "code 1",
  "06.03.2023", "code 1",
  "07.03.2023", "code 1",
  "08.03.2023", "code 1",
  "09.03.2023", "code 0.3",
  "10.03.2023", "code 0.8",
  "11.03.2023", "code 0.4",
  "10.04.2023", "code 0.7",
  "11.04.2023", "code 0.3",
  "14.04.2023", "code 0.6",
  "15.04.2023", "code 1",
  "16.04.2023", "code 1",
  "17.04.2023", "code 1",
  "18.04.2023", "code 1",
  "19.04.2023", "code 1",
  "20.04.2023", "code 1",
  "21.04.2023", "code 1",
  "22.04.2023", "code 1",
  "23.04.2023", "code 1",
  "24.04.2023", "code 1",
  "25.04.2023", "code 1",
  "26.04.2023", "code 0.6",
  "27.04.2023", "code 0.4",
  "11.05.2023", "code 1",
  "12.05.2023", "code 1",
  "13.05.2023", "code 1",
  "16.05.2023", "code 1",
  "17.05.2023", "code 1",
  "18.05.2023", "code 1",
  "19.05.2023", "code 1",
  "20.05.2023", "code 1",
  "21.05.2023", "code 0.6",
  "22.05.2023", "code 0.4",
  "29.05.2023", "code 1",
  "30.05.2023", "code 1",
  "31.05.2023", "code 1",
  "09.06.2023", "code 1",
  "10.06.2023", "code 1",
  "11.06.2023", "code 1",
  "12.06.2023", "code 1",
  "19.06.2023", "code 0.6",
  "20.06.2023", "code 1",
  "21.06.2023", "code 1",
  "22.06.2023", "code 1",
  "23.06.2023", "code 1",
  "24.06.2023", "code 1",
  "25.06.2023", "code 1",
  "26.06.2023", "code 1",
  "27.06.2023", "code 0.6",
  "28.06.2023", "code 1",
  "29.06.2023", "code 1",
  "30.06.2023", "code 1",
  "09.07.2023", "code 1",
  "10.07.2023", "code 1",
  "11.07.2023", "code 1",
  "12.07.2023", "code 1",
  "13.07.2023", "code 1",
  "14.07.2023", "code 1",
  "20.07.2023", "code 0.6",
  "21.07.2023", "code 1",
  "22.07.2023", "code 1",
  "23.07.2023", "code 1",
  "24.07.2023", "code 1",
  "25.07.2023", "code 1",
  "26.07.2023", "code 1",
  "27.07.2023", "code 1",
  "28.07.2023", "code 1",
  "29.07.2023", "code 1",
  "30.07.2023", "code 1",
  "31.07.2023", "code 1",
  "01.08.2023", "code 1",
  "16.08.2023", "code 1",
  "17.08.2023", "code 1",
  "18.08.2023", "code 1",
  "19.08.2023", "code 1",
  "22.08.2023", "code 1",
  "23.08.2023", "code 1",
  "24.08.2023", "code 1",
  "25.08.2023", "code 1",
  "26.08.2023", "code 1",
  "27.08.2023", "code 1",
  "28.08.2023", "code 1",
  "29.08.2023", "code 1",
  "30.08.2023", "code 1",
  "31.08.2023", "code 1",
  "01.09.2023", "code 1",
  "02.09.2023", "code 1",
  "03.09.2023", "code 1",
  "04.09.2023", "code 1",
  "05.09.2023", "code 1",
  "06.09.2023", "code 1",
  "07.09.2023", "code 1",
  "08.09.2023", "code 1",
  "16.10.2023", "code 0.4"
) %>%
  mutate(date = lubridate::dmy(date_chr),
         code = readr::parse_number(detail)) %>%
  dplyr::select(date, code)


# JC1 2024 grazing detail: date -> code (+ optional explicit HH:MM times)
jc1_grazing_detail_2024 <- tibble::tribble(
  ~date_chr,    ~detail,
  "27.02.2024", "code 0.5",
  "29.02.2024", "code 0.2",
  "03.03.2024", "code 0.5",
  "04.03.2024", "code 0.2",
  "05.03.2024", "code 0.2",   # explicit morning + also treat as 0.5 for evening
  "05.03.2024", "code 0.5",
  "06.03.2024", "code 0.5",
  "07.03.2024", "code 0.5",
  "08.03.2024", "code 0.5",
  "09.03.2024", "code 0.5",
  "11.03.2024", "code 0.5",
  "12.03.2024", "code 0.2",
  "16.04.2024", "code 0.5",
  "17.04.2024", "code 0.5",
  "18.04.2024", "code 1",
  "19.04.2024", "code 1",
  "20.04.2024", "code 1",
  "21.04.2024", "code 1",
  "22.04.2024", "code 1",
  "23.04.2024", "code 1",
  "24.04.2024", "code 1",
  "25.04.2024", "code 1",
  "26.04.2024", "code 1",
  "27.04.2024", "code 1",
  "28.04.2024", "code 1",
  "29.04.2024", "code 1",
  "30.04.2024", "code 0.8",
  "01.05.2024", "code 0.6",
  "06.05.2024", "code 0.6",
  "07.05.2024", "code 1",
  "08.05.2024", "code 1",
  "10.05.2024", "code 1",
  "17.05.2024", "code 1",
  "18.05.2024", "code 1",
  "19.05.2024", "code 1",
  "20.05.2024", "code 1",
  "21.05.2024", "code 1",
  "22.05.2024", "code 1",
  "23.05.2024", "code 1",
  "24.05.2024", "code 1",
  "25.05.2024", "code 0.4",
  "26.05.2024", "code 0.6",
  "27.05.2024", "code 1",
  "28.05.2024", "code 0.4",
  "09.06.2024", "code 0.6",
  "10.06.2024", "code 1",
  "11.06.2024", "code 1",
  "12.06.2024", "code 1",
  "13.06.2024", "code 1",
  "14.06.2024", "code 1",
  "15.06.2024", "code 1",
  "16.06.2024", "code 1",
  "17.06.2024", "code 1",
  "18.06.2024", "code 1",
  "19.06.2024", "code 1",
  "20.06.2024", "code 1",
  "21.06.2024", "code 1",
  "22.06.2024", "code 1",
  "23.06.2024", "code 1",
  "24.06.2024", "code 1",
  "25.06.2024", "code 1",
  "26.06.2024", "code 1",
  "27.06.2024", "code 1",
  "23.07.2024", "code 1",
  "24.07.2024", "code 1",
  "25.07.2024", "code 1",
  "26.07.2024", "code 1",
  "27.07.2024", "code 1",
  "28.07.2024", "code 1",
  "29.07.2024", "code 1",
  "30.07.2024", "code 1",
  "31.07.2024", "code 1",
  "01.08.2024", "code 1",
  "02.08.2024", "code 1",
  "03.08.2024", "code 1",
  "04.08.2024", "code 1",
  "05.08.2024", "code 1",
  "06.08.2024", "code 1",
  "07.08.2024", "code 1",
  "08.08.2024", "code 1",
  "09.08.2024", "code 1",
  "25.09.2024", "code 1",
  "26.09.2024", "code 0.4",
  "27.09.2024", "code 1",
  "28.09.2024", "code 0.4",
  "12.10.2024", "code 1",
  "13.10.2024", "code 1"
) %>%
  mutate(
    date = lubridate::dmy(date_chr),
    date = dplyr::case_when(
      lubridate::year(date) %in% c(2025, 2026, 2027) ~
        lubridate::make_date(2024, lubridate::month(date), lubridate::day(date)),
      TRUE ~ date
    ),
    code  = readr::parse_number(detail),
    times = stringr::str_extract(detail, "\\b\\d{1,2}:\\d{2}-\\d{1,2}:\\d{2}\\b"),
    times = stringr::str_replace_all(times, "^(\\d):(\\d{2})", "0\\1:\\2"),
    times = stringr::str_replace_all(times, "-(\\d):(\\d{2})$", "-0\\1:\\2"),
    t_start = dplyr::if_else(!is.na(times), stringr::str_sub(times, 1, 5), NA_character_),
    t_end   = dplyr::if_else(!is.na(times), stringr::str_sub(times, 7, 11), NA_character_)
  ) %>%
  dplyr::select(date, code, t_start, t_end)

# JC1 2020 grazing detail: date -> grazing code
jc1_grazing_detail_2020 <- tibble::tribble(
  ~date_chr,    ~detail,
  "04.02.2020", "code 0.3",
  "05.02.2020", "code 0.7",
  "06.02.2020", "code 0.7",
  "07.02.2020", "code 0.3",
  "10.02.2020", "code 0.2",
  "03.03.2020", "code 0.8",
  "04.03.2020", "code 0.3",
  "05.03.2020", "code 0.8",
  "06.03.2020", "code 0.8",
  "07.03.2020", "code 0.4",
  "09.03.2020", "code 0.2",
  "12.03.2020", "code 0.3",
  "13.03.2020", "code 0.4",
  "16.03.2020", "code 0.8",
  "17.03.2020", "code 0.4",
  "19.03.2020", "code 0.8",
  "20.03.2020", "code 1",
  "21.03.2020", "code 1",
  "22.03.2020", "code 0.4",
  "10.04.2020", "code 0.6",
  "11.04.2020", "code 1",
  "12.04.2020", "code 1",
  "13.04.2020", "code 1",
  "14.04.2020", "code 1",
  "15.04.2020", "code 1",
  "16.04.2020", "code 1",
  "17.04.2020", "code 1",
  "18.04.2020", "code 0.4",
  "03.05.2020", "code 0.6",
  "04.05.2020", "code 1",
  "05.05.2020", "code 1",
  "06.05.2020", "code 0.6",
  "06.05.2020", "code 0.4",
  "07.05.2020", "code 1",
  "08.05.2020", "code 1",
  "09.05.2020", "code 0.4",
  "25.05.2020", "code 1",
  "26.05.2020", "code 1",
  "27.05.2020", "code 1",
  "28.05.2020", "code 1",
  "29.05.2020", "code 0.6",
  "29.05.2020", "code 0.4",
  "30.05.2020", "code 1",
  "31.05.2020", "code 0.6",
  "31.05.2020", "code 0.4",
  "01.06.2020", "code 1",
  "02.06.2020", "code 1",
  "17.06.2020", "code 0.6",
  "18.06.2020", "code 1",
  "19.06.2020", "code 1",
  "20.06.2020", "code 1",
  "21.06.2020", "code 0.6",
  "21.06.2020", "code 0.4",
  "22.06.2020", "code 1",
  "23.06.2020", "code 1",
  "24.06.2020", "code 0.4",
  "25.06.2020", "code 0.4",
  "09.07.2020", "code 1",
  "10.07.2020", "code 1",
  "11.07.2020", "code 1",
  "12.07.2020", "code 1",
  "13.07.2020", "code 1",
  "14.07.2020", "code 1",
  "15.07.2020", "code 0.6",
  "15.07.2020", "code 0.4",
  "16.07.2020", "code 1",
  "17.07.2020", "code 1",
  "18.07.2020", "code 0.4",
  "01.08.2020", "code 1",
  "02.08.2020", "code 1",
  "03.08.2020", "code 1",
  "04.08.2020", "code 1",
  "05.08.2020", "code 1",
  "06.08.2020", "code 1",
  "07.08.2020", "code 1",
  "08.08.2020", "code 1",
  "09.08.2020", "code 1",
  "10.08.2020", "code 1",
  "11.08.2020", "code 1",
  "12.08.2020", "code 0.4",
  "01.09.2020", "code 1",
  "02.09.2020", "code 1",
  "03.09.2020", "code 1",
  "04.09.2020", "code 1",
  "05.09.2020", "code 1",
  "06.09.2020", "code 1",
  "07.09.2020", "code 1",
  "08.09.2020", "code 1",
  "09.09.2020", "code 1",
  "10.09.2020", "code 1",
  "11.09.2020", "code 1",
  "12.09.2020", "code 1",
  "13.09.2020", "code 1",
  "14.09.2020", "code 1",
  "15.09.2020", "code 1",
  "16.09.2020", "code 1",
  "17.09.2020", "code 1",
  "18.09.2020", "code 1",
  "19.09.2020", "code 1",
  "20.09.2020", "code 1",
  "21.09.2020", "code 0.4",
  "21.10.2020", "code 0.8",
  "22.10.2020", "code 0.8",
  "23.10.2020", "code 1",
  "24.10.2020", "code 1",
  "25.10.2020", "code 0.8",
  "26.10.2020", "code 0.8",
  "27.10.2020", "code 0.4",
  "28.10.2020", "code 0.8",
  "29.10.2020", "code 0.3",
  "30.10.2020", "code 0.8",
  "02.11.2020", "code 0.4",
  "03.11.2020", "code 0.8",
  "04.11.2020", "code 0.8",
  "05.11.2020", "code 0.8",
  "06.11.2020", "code 0.8",
  "07.11.2020", "code 0.4"
) %>%
  mutate(date = lubridate::dmy(date_chr),
         code = readr::parse_number(detail)) %>%
  dplyr::select(date, code)

# JC2 2023 grazing detail: date -> grazing code
jc2_grazing_detail_2023 <- tibble::tribble(
  ~date_chr,    ~detail,
  "21.03.2023", "code 0.3",
  "22.03.2023", "code 0.3",
  "23.03.2023", "code 0.3",
  "24.03.2023", "code 0.8",
  "24.04.2023", "code 1",
  "25.04.2023", "code 1",
  "21.05.2023", "code 0.6",
  "22.05.2023", "code 1",
  "23.05.2023", "code 1",
  "24.05.2023", "code 1",
  "21.06.2023", "code 0.6",
  "22.06.2023", "code 1",
  "23.06.2023", "code 1",
  "24.06.2023", "code 1",
  "25.06.2023", "code 0.4",
  "21.07.2023", "code 1",
  "22.07.2023", "code 1",
  "17.08.2023", "code 1",
  "18.08.2023", "code 1",
  "19.08.2023", "code 1",
  "20.08.2023", "code 1",
  "21.08.2023", "code 0.4",
  "23.09.2023", "code 0.4",
  "24.09.2023", "code 0.4",
  "25.09.2023", "code 0.4",
  "26.09.2023", "code 0.4",
  "27.09.2023", "code 0.4",
  "28.09.2023", "code 0.8",
  "29.09.2023", "code 0.4",
  "03.10.2023", "code 0.4",
  "04.10.2023", "code 0.4"
) %>%
  mutate(date = lubridate::dmy(date_chr),
         code = readr::parse_number(detail)) %>%
  dplyr::select(date, code)


# JC2 2024 grazing detail: date -> code (+ optional explicit HH:MM times)
jc2_grazing_detail_2024 <- tibble::tribble(
  ~date_chr,    ~detail,
  "07.02.2024", "11:00-14:00",
  "12.02.2024", "12:00-15:00",
  "13.02.2024", "09:00-12:00",
  "15.04.2024", "code 0.2",
  "16.04.2024", "code 0.2",
  "17.04.2024", "code 0.2",
  "18.04.2024", "code 0.3",
  "22.05.2024", "code 0.6",
  "23.05.2024", "code 0.6",
  "24.05.2024", "code 0.4",
  "24.05.2024", "code 0.6",
  "25.05.2024", "code 0.6",
  "27.05.2024", "code 1",
  "27.05.2024", "code 1",
  "28.05.2024", "code 1",
  "20.06.2024", "code 1",
  "21.06.2024", "code 1",
  "22.06.2024", "code 1",
  "23.06.2024", "code 1",
  "24.06.2024", "code 1",
  "25.06.2024", "code 0.6",
  "26.07.2024", "code 1",
  "27.07.2024", "code 1",
  "27.07.2024", "code 1",
  "28.07.2024", "code 1",
  "28.07.2024", "code 1",
  "29.07.2024", "code 1",
  "30.07.2024", "code 1",
  "09.09.2024", "code 0.4",
  "10.09.2024", "code 0.4",
  "11.09.2024", "code 0.4",
  "12.09.2024", "code 0.4",
  "13.09.2024", "code 1",
  "14.09.2024", "code 1",
  "15.09.2024", "code 1"
) %>%
  mutate(
    date = lubridate::dmy(date_chr),
    date = dplyr::case_when(
      lubridate::month(date) == 9 & lubridate::year(date) > 2024 ~
        lubridate::make_date(2024, 9, lubridate::day(date)),
      TRUE ~ date
    ),
    code  = readr::parse_number(detail),
    times = stringr::str_extract(detail, "\\b\\d{2}:\\d{2}-\\d{2}:\\d{2}\\b"),
    t_start = dplyr::if_else(!is.na(times), stringr::str_sub(times, 1, 5), NA_character_),
    t_end   = dplyr::if_else(!is.na(times), stringr::str_sub(times, 7, 11), NA_character_)
  ) %>%
  dplyr::select(date, code, t_start, t_end)

# de-duplicate the 2024 grazing date lists (they contain repeats)
management[["JC1"]]$events[["2024"]][["Grazing"]] <- unique(as.Date(management[["JC1"]]$events[["2024"]][["Grazing"]]))
management[["JC2"]]$events[["2024"]][["Grazing"]] <- unique(as.Date(management[["JC2"]]$events[["2024"]][["Grazing"]]))

# pad start/end times to HH:MM (24:00 kept as-is, no 23:59)
pad_times <- function(df) {
  df %>%
    mutate(
      start_time = stringr::str_pad(start_time, 5, pad = "0"),
      end_time   = stringr::str_pad(end_time,   5, pad = "0")
    )
}

# expand code rows into 1-3 daily windows; same-day windows share one occ_group
expand_code_windows_grouped <- function(df_code, allowed_next_day = NULL, event_label = "Grazing") {
  dfc <- df_code %>%
    dplyr::left_join(code_map, by = "code") %>%
    dplyr::mutate(origin_date = date,
                  occ_group   = paste0("occ_", format(origin_date, "%Y-%m-%d")))
  
  w <- dplyr::bind_rows(
    dfc %>% dplyr::filter(!is.na(t1_start)) %>%
      dplyr::transmute(date, event = event_label, start_time = t1_start, end_time = t1_end, occ_group),
    dfc %>% dplyr::filter(!is.na(t2_start)) %>%
      dplyr::transmute(date, event = event_label, start_time = t2_start, end_time = t2_end, occ_group),
    dfc %>% dplyr::filter(!is.na(t3_start)) %>%
      dplyr::transmute(date, event = event_label, start_time = t3_start, end_time = t3_end, occ_group)
  )
  
  # 0.6 carries into next day 00:00–07:30 (same occ_group!)
  carry <- dfc %>%
    dplyr::filter(carry_next) %>%
    dplyr::transmute(
      date       = origin_date + lubridate::days(1),
      event      = event_label,
      start_time = "00:00",
      end_time   = "07:30",
      occ_group
    )
  
  if (!is.null(allowed_next_day)) {
    carry <- carry %>% dplyr::filter(date %in% allowed_next_day)
  }
  
  dplyr::bind_rows(w, carry) %>% pad_times()
}

# number occurrences chronologically within each event category (Grazing_1, ...)
number_events_all <- function(ev_tbl, tz = "UTC") {
  if (!nrow(ev_tbl)) return(ev_tbl)
  
  # Normalise & prepare
  ev_tbl <- ev_tbl %>%
    pad_times() %>%
    mutate(
      event    = as.character(event),
      category = stringr::str_remove(event, "_\\d+$"),   # base label (e.g., "Grazing")
      has_idx  = stringr::str_detect(event, "_\\d+$")    # already numbered?
    )
  
  # Ensure an occ_group exists for grouping occurrences:
  # - keep existing occ_group if supplied (e.g., from code-expansion for Grazing)
  # - for already-numbered rows, use a stable key so they keep their label
  # - otherwise: one occurrence per (category, day)
  if (!"occ_group" %in% names(ev_tbl)) {
    ev_tbl$occ_group <- NA_character_
  }
  ev_tbl <- ev_tbl %>%
    mutate(
      occ_group = dplyr::case_when(
        has_idx ~ paste0(category, "::", event),  # keep pre-numbered labels intact
        is.na(occ_group) | occ_group == "" ~ paste0(category, "::", format(as.Date(date), "%Y-%m-%d")),
        TRUE ~ occ_group
      ),
      start_ct = as.POSIXct(paste(date, start_time), tz = tz)
    )
  
  # Order groups in time within each category and assign indices
  grp_order <- ev_tbl %>%
    group_by(category, occ_group) %>%
    summarise(group_start = min(start_ct, na.rm = TRUE), .groups = "drop") %>%
    arrange(category, group_start) %>%
    group_by(category) %>%
    mutate(idx = dplyr::row_number()) %>%
    ungroup()
  
  # Build final numbered labels
  ev_tbl_num <- ev_tbl %>%
    left_join(grp_order, by = c("category","occ_group")) %>%
    mutate(
      event = if_else(has_idx, event, paste0(category, "_", idx))
    ) %>%
    dplyr::select(date, event, start_time, end_time) %>%
    arrange(as.Date(date), event, start_time)
  
  ev_tbl_num
}


# JC1 2020: grazing (code windows) + fixed-window fertiliser events
allowed_jc1_2020 <- as.Date(management[["JC1"]]$events[["2020"]][["Grazing"]])

jc1_grazing_windows_2020 <- expand_code_windows_grouped(
  jc1_grazing_detail_2020,
  allowed_next_day = allowed_jc1_2020
)

jc1_fertiliser_windows_2020 <- tibble::tibble(
  date       = as.Date(management[["JC1"]]$events[["2020"]][["Fertiliser"]]),
  event      = "Fertiliser",
  start_time = "10:30",
  end_time   = "11:00"
) %>% pad_times()

jc1_events_2020 <- dplyr::bind_rows(
  jc1_grazing_windows_2020,
  jc1_fertiliser_windows_2020
) %>% number_events_all()

# JC1 2023: grazing (code windows) + fixed-window fertiliser events
allowed_jc1_2023 <- as.Date(management[["JC1"]]$events[["2023"]][["Grazing"]])

jc1_grazing_windows_2023 <- expand_code_windows_grouped(
  jc1_grazing_detail_2023,
  allowed_next_day = allowed_jc1_2023
)

jc1_fertiliser_windows_2023 <- tibble::tibble(
  date       = as.Date(management[["JC1"]]$events[["2023"]][["Fertiliser"]]),
  event      = "Fertiliser",
  start_time = "10:30",
  end_time   = "11:00"
) %>% pad_times()

jc1_events_2023 <- dplyr::bind_rows(
  jc1_grazing_windows_2023,
  jc1_fertiliser_windows_2023
) %>% number_events_all()
# JC1 2024: explicit-time rows kept per row, code rows grouped and expanded
# explicit rows keep their own group per row
gw_explicit_jc1_2024 <- jc1_grazing_detail_2024 %>%
  dplyr::filter(!is.na(t_start)) %>%
  dplyr::transmute(
    date, event = "Grazing",
    start_time = t_start, end_time = t_end,
    occ_group  = paste0("occ_", format(date, "%Y-%m-%d"))  # <-- same pattern
  ) %>% pad_times()


# Code rows use grouped expansion (one occ_group per origin date)
gw_code_base_jc1_2024 <- jc1_grazing_detail_2024 %>%
  dplyr::filter(is.na(t_start) & !is.na(code)) %>%
  dplyr::select(date, code)

allowed_jc1_2024 <- as.Date(management[["JC1"]]$events[["2024"]][["Grazing"]])
gw_code_jc1_2024 <- expand_code_windows_grouped(gw_code_base_jc1_2024, allowed_next_day = allowed_jc1_2024)

jc1_fertiliser_windows_2024 <- tibble::tibble(
  date       = as.Date(management[["JC1"]]$events[["2024"]][["Fertiliser"]]),
  event      = "Fertiliser",
  start_time = "10:30",
  end_time   = "11:00"
) %>% pad_times()

jc1_events_2024 <- dplyr::bind_rows(
  gw_explicit_jc1_2024,
  gw_code_jc1_2024,
  jc1_fertiliser_windows_2024
) %>% number_events_all()

# JC2 2023: grazing (code windows) + fixed-window fertiliser events
jc2_grazing_detail_2023 <- jc2_grazing_detail_2023 %>% dplyr::select(date, code)
allowed_jc2_2023 <- as.Date(management[["JC2"]]$events[["2023"]][["Grazing"]])

jc2_grazing_windows_2023 <- expand_code_windows_grouped(jc2_grazing_detail_2023, allowed_next_day = allowed_jc2_2023)

jc2_fertiliser_windows_2023 <- tibble::tibble(
  date  = as.Date(management[["JC2"]]$events[["2023"]][["Fertiliser"]]),
  event = "Fertiliser",
  start_time = "10:30",
  end_time   = "11:00"
) %>% pad_times()

jc2_events_2023 <- dplyr::bind_rows(
  jc2_grazing_windows_2023,
  jc2_fertiliser_windows_2023
) %>% number_events_all()

# JC2 2024: explicit-time rows kept per row, code rows grouped and expanded
gw_explicit_jc2_2024 <- jc2_grazing_detail_2024 %>%
  dplyr::filter(!is.na(t_start)) %>%
  dplyr::transmute(
    date, event = "Grazing",
    start_time = t_start, end_time = t_end,
    occ_group  = paste0("occ_", format(date, "%Y-%m-%d"))  # <-- same pattern
  ) %>% pad_times()


gw_code_base_jc2_2024 <- jc2_grazing_detail_2024 %>%
  dplyr::filter(is.na(t_start) & !is.na(code)) %>%
  dplyr::select(date, code)

allowed_jc2_2024 <- as.Date(management[["JC2"]]$events[["2024"]][["Grazing"]])
gw_code_jc2_2024 <- expand_code_windows_grouped(gw_code_base_jc2_2024, allowed_next_day = allowed_jc2_2024)

jc2_grazing_windows_2024 <- dplyr::bind_rows(gw_explicit_jc2_2024, gw_code_jc2_2024)

# Fertiliser & slurry
jc2_fert_dates_raw      <- as.Date(management[["JC2"]]$events[["2024"]][["Fertiliser"]])
jc2_fertiliser_dates_24 <- setdiff(jc2_fert_dates_raw, as.Date("2024-04-19"))
jc2_slurry_dates_24     <- as.Date("2024-04-19")

jc2_fertiliser_windows_2024 <- tibble::tibble(
  date  = jc2_fertiliser_dates_24,
  event = "Fertiliser",
  start_time = "10:30",
  end_time   = "11:00"
) %>% pad_times()

jc2_slurry_windows_2024 <- tibble::tibble(
  date  = jc2_slurry_dates_24,
  event = "Slurry",
  start_time = "10:30",
  end_time   = "11:00"
) %>% pad_times()

jc2_events_2024 <- dplyr::bind_rows(
  jc2_grazing_windows_2024,
  jc2_fertiliser_windows_2024,
  jc2_slurry_windows_2024
) %>% number_events_all()

# sort each event table and drop exact duplicate rows
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

# finalize all JC1 / JC2 site-year event tables
jc1_events_2020 <- finalize_site_events(jc1_events_2020)
jc1_events_2023 <- finalize_site_events(jc1_events_2023)
jc1_events_2024 <- finalize_site_events(jc1_events_2024)

jc2_events_2023 <- finalize_site_events(jc2_events_2023)
jc2_events_2024 <- finalize_site_events(jc2_events_2024)
