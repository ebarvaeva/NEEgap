# monthly_temperature_rainfall_table.R — Monthly rain & temperature summary from Met Eireann hourly data
#
# Reads the hourly met CSV, keeps 2020 / 2023 / 2024, sums rainfall and averages
# temperature per calendar month, then reshapes to one row per month with a
# year x variable column for each value. Prints the 12-row table; nothing saved.
#
# Input:  data/meteireann_data/met_hourly.csv (columns: date "dmy hm", rain, temp)
# Output: none — the wide monthly table is printed to the console.

library(readr); library(dplyr); library(lubridate); library(tidyr)

CSV <- "data/meteireann_data/met_hourly.csv"

# Parse the "dd/mm/yyyy hh:mm" date, coerce numerics, keep the three study years.
df <- read_csv(CSV, na = c("", " ", "NA"), show_col_types = FALSE) |>
  mutate(
    datetime = dmy_hm(date),
    rain     = as.numeric(rain),
    temp     = as.numeric(temp),
    year     = year(datetime),
    month    = month(datetime)
  ) |>
  filter(year %in% c(2020, 2023, 2024))

# Monthly rainfall total and mean temperature.
monthly <- df |>
  group_by(year, month) |>
  summarise(
    rain_mm = sum(rain, na.rm = TRUE),
    temp_C  = mean(temp, na.rm = TRUE),
    .groups = "drop"
  )

# Wide format: one row per month, one column per year x variable, rounded to 1 dp.
wide <- monthly |>
  pivot_wider(
    names_from  = year,
    values_from = c(temp_C, rain_mm),
    names_glue  = "{.value}_{year}"
  ) |>
  mutate(month_name = month.name[month]) |>
  select(month_name,
         temp_C_2020, temp_C_2023, temp_C_2024,
         rain_mm_2020, rain_mm_2023, rain_mm_2024) |>
  mutate(across(where(is.numeric), \(x) round(x, 1)))

print(wide, n = 12)
