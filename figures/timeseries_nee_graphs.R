# timeseries_nee.R — JC1-vs-JC2 NEE figures for the paper (Figs 1, 2, 3)
#
# Standalone plotting script. Loads both sites' prepared NEE series and draws:
#   Fig 1 — NEE time series per site with grazing/fertiliser event lines
#           (JC1 2020, JC1 2023-2024, JC2 2023-2024 stacked)
#   Fig 2 — JC1-vs-JC2 NEE scatter, 2x2 facets by JC1 grazing window,
#           coloured by JC1 days-since-fertiliser
#   Fig 3 — same as Fig 2 but faceted/coloured by JC2 management
#   Fig 4 — Figs 2 and 3 side by side
#
# Inputs (in data/data_prepared/):
#   JC1.rds, JC2.rds
#   required columns: timestamp, NEE_orig, Grazing_days_since,
#                     Fertiliser_days_since
#
# Outputs (in graphs/nee_sites/):
#   fig1_nee_timeseries.{png,pdf}
#   fig2_nee_jc1_vs_jc2_grazing_windows_jc1mgmt.{png,pdf}
#   fig3_nee_jc1_vs_jc2_grazing_windows_jc2mgmt.{png,pdf}
#   fig4_scatter_combined.png
#
# To adapt: point DATA_DIR at the folder holding {SITE}.rds and adjust SITES.

DATA_DIR <- here::here("data", "data_prepared")
OUT_DIR  <- here::here("graphs", "nee_sites")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SITES      <- c("JC1", "JC2")
SITE_COLS  <- c(JC1 = "#377EB8", JC2 = "#E41A1C")   # blue / red

# Management window: rows within N days of an event are "post-event"
POST_GRAZ_DAYS <- 30
POST_FERT_DAYS <- 14

# Season order (Irish meteorological)
SEASON_LEVELS <- c("Winter", "Spring", "Summer", "Autumn")

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(lubridate); library(scales); library(patchwork)
})

# save_fig(): write one plot as both PNG and PDF.
save_fig <- function(p, name, w = 14, h = 7) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p,
         width = w, height = h)
}

base_theme <- theme_bw(base_size = 16) +
  theme(
    panel.grid.minor  = element_blank(),
    legend.position   = "bottom",
    plot.background   = element_rect(fill = "white", colour = NA),
    strip.background  = element_rect(fill = "grey92"),
    strip.text        = element_text(face = "bold", size = 15),
    axis.text         = element_text(size = 14),
    axis.title        = element_text(size = 15),
    legend.text       = element_text(size = 13),
    legend.title      = element_text(size = 14)
  )

# Load both sites, add time parts and grazing/fertiliser phase labels.
df_all <- map_dfr(SITES, function(s) {
  path <- file.path(DATA_DIR, paste0(s, ".rds"))
  if (!file.exists(path)) { return(NULL) }
  readRDS(path) %>% mutate(site = s)
}) %>%
  mutate(
    timestamp = as.POSIXct(timestamp, tz = "UTC"),
    site      = factor(site, levels = SITES),
    month     = month(timestamp),
    year      = year(timestamp),
    doy       = yday(timestamp),
    season    = case_when(
      month %in% c(11, 12, 1)  ~ "Winter",
      month %in% c(2, 3, 4)    ~ "Spring",
      month %in% c(5, 6, 7)    ~ "Summer",
      month %in% c(8, 9, 10)   ~ "Autumn"
    ),
    season = factor(season, levels = SEASON_LEVELS),
    # Management phase categories
    graz_phase = case_when(
      !is.na(Grazing_days_since) & Grazing_days_since == 0              ~ "During grazing",
      !is.na(Grazing_days_since) & Grazing_days_since <= POST_GRAZ_DAYS ~ "Post-grazing (\u226430 d)",
      !is.na(Grazing_days_since) & Grazing_days_since >  POST_GRAZ_DAYS ~ "Recovery (>30 d)",
      TRUE ~ NA_character_
    ),
    graz_phase = factor(graz_phase,
                        levels = c("During grazing",
                                   "Post-grazing (\u226430 d)",
                                   "Recovery (>30 d)")),
    fert_phase = case_when(
      !is.na(Fertiliser_days_since) & Fertiliser_days_since == 0               ~ "During fertilisation",
      !is.na(Fertiliser_days_since) & Fertiliser_days_since <= POST_FERT_DAYS  ~ "Post-fertilisation (\u226414 d)",
      !is.na(Fertiliser_days_since) & Fertiliser_days_since >  POST_FERT_DAYS  ~ "Background (>14 d)",
      TRUE ~ NA_character_
    ),
    fert_phase = factor(fert_phase,
                        levels = c("During fertilisation",
                                   "Post-fertilisation (\u226414 d)",
                                   "Background (>14 d)"))
  )


# --- FIGURE 1 — NEE time series with grazing (green) / fertiliser (brown) lines ---

# JC1 2020: subset and pull grazing/fertiliser event timestamps.
d1 <- df_all %>%
  filter(site == "JC1",
         year(timestamp) == 2020)

ts_graz1 <- d1 %>%
  filter(!is.na(Grazing_days_since),
         Grazing_days_since == 0) %>%
  pull(timestamp)

ts_fert1 <- d1 %>%
  filter(!is.na(Fertiliser_days_since),
         Fertiliser_days_since == 0) %>%
  pull(timestamp)

p1a <- ggplot(
  d1 %>% filter(!is.na(NEE_orig)),
  aes(timestamp, NEE_orig)
) +
  geom_vline(
    xintercept = as.numeric(ts_graz1),
    colour = "#4DAF4A",
    linewidth = 0.35,
    linetype = "dashed",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = as.numeric(ts_fert1),
    colour = "#A65628",
    linewidth = 0.85,
    linetype = "dashed"
  ) +
  geom_point(colour = "grey20", size = 0.3, alpha = 0.4) +
  scale_x_datetime(
    date_breaks = "1 month",
    labels = function(x) {
      ifelse(
        format(x, "%m") == "01",
        paste0(format(x, "%b"), "\n", format(x, "%Y")),
        format(x, "%b")
      )
    }
  ) +
  scale_y_continuous(
    limits = c(-40, 40),
    name = expression(NEE~(mu*mol~m^{-2}~s^{-1}))
  ) +
  labs(title = "JC1 2020", x = NULL) +
  base_theme


# JC1 2023-2024: same construction.
d2 <- df_all %>%
  filter(site == "JC1",
         year(timestamp) %in% c(2023, 2024))

ts_graz2 <- d2 %>%
  filter(!is.na(Grazing_days_since),
         Grazing_days_since == 0) %>%
  pull(timestamp)

ts_fert2 <- d2 %>%
  filter(!is.na(Fertiliser_days_since),
         Fertiliser_days_since == 0) %>%
  pull(timestamp)

p1b <- ggplot(
  d2 %>% filter(!is.na(NEE_orig)),
  aes(timestamp, NEE_orig)
) +
  geom_vline(
    xintercept = as.numeric(ts_graz2),
    colour = "#4DAF4A",
    linewidth = 0.35,
    linetype = "dashed",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = as.numeric(ts_fert2),
    colour = "#A65628",
    linewidth = 0.85,
    linetype = "dashed"
  ) +
  geom_point(colour = "grey20", size = 0.3, alpha = 0.4) +
  scale_x_datetime(
    date_breaks = "1 month",
    labels = function(x) {
      ifelse(
        format(x, "%m") == "01",
        paste0(format(x, "%b"), "\n", format(x, "%Y")),
        format(x, "%b")
      )
    }
  ) +
  scale_y_continuous(
    limits = c(-40, 40),
    name = expression(NEE~(mu*mol~m^{-2}~s^{-1}))
  ) +
  labs(title = "JC1 2023–2024", x = NULL) +
  base_theme


# JC2 2023-2024: same construction.
d3 <- df_all %>%
  filter(site == "JC2",
         year(timestamp) %in% c(2023, 2024))

ts_graz3 <- d3 %>%
  filter(!is.na(Grazing_days_since),
         Grazing_days_since == 0) %>%
  pull(timestamp)

ts_fert3 <- d3 %>%
  filter(!is.na(Fertiliser_days_since),
         Fertiliser_days_since == 0) %>%
  pull(timestamp)

p1c <- ggplot(
  d3 %>% filter(!is.na(NEE_orig)),
  aes(timestamp, NEE_orig)
) +
  geom_vline(
    xintercept = as.numeric(ts_graz3),
    colour = "#4DAF4A",
    linewidth = 0.35,
    linetype = "dashed",
    alpha = 0.8
  ) +
  geom_vline(
    xintercept = as.numeric(ts_fert3),
    colour = "#A65628",
    linewidth = 0.85,
    linetype = "dashed"
  ) +
  geom_point(colour = "grey20", size = 0.3, alpha = 0.4) +
  scale_x_datetime(
    date_breaks = "1 month",
    labels = function(x) {
      ifelse(
        format(x, "%m") == "01",
        paste0(format(x, "%b"), "\n", format(x, "%Y")),
        format(x, "%b")
      )
    }
  ) +
  scale_y_continuous(
    limits = c(-40, 40),
    name = expression(NEE~(mu*mol~m^{-2}~s^{-1}))
  ) +
  labs(title = "JC2 2023–2024", x = NULL) +
  base_theme


# Stack the three panels and save Fig 1.
p1 <- p1a / p1b / p1c

save_fig(
  p1,
  "fig1_nee_timeseries",
  w = 15,
  h = 11
)


# --- FIGURE 2 — JC1-vs-JC2 NEE scatter, faceted by JC1 grazing window, coloured by JC1 dsf ---

# Wide table: one row per shared timestamp, NEE + management for both sites.
df_wide_base <- df_all %>%
  filter(!is.na(NEE_orig)) %>%
  select(timestamp, site, NEE_orig, Grazing_days_since, Fertiliser_days_since) %>%
  pivot_wider(
    id_cols     = timestamp,
    names_from  = site,
    values_from = c(NEE_orig, Grazing_days_since, Fertiliser_days_since)
  ) %>%
  filter(!is.na(NEE_orig_JC1), !is.na(NEE_orig_JC2))

# Assign each row to a JC1 grazing-recovery window.
df_wide_jc1 <- df_wide_base %>%
  filter(!is.na(Grazing_days_since_JC1),
         !is.na(Fertiliser_days_since_JC1)) %>%
  mutate(
    graz_window = case_when(
      Grazing_days_since_JC1 <= 1                                  ~ "0\u20131 d",
      Grazing_days_since_JC1 >  1 & Grazing_days_since_JC1 <= 7   ~ "1\u20137 d",
      Grazing_days_since_JC1 >  7 & Grazing_days_since_JC1 <= 30  ~ "7\u201330 d",
      Grazing_days_since_JC1 >  30                                 ~ "> 30 d",
      TRUE ~ NA_character_
    ),
    graz_window = factor(graz_window,
                         levels = c("0\u20131 d", "1\u20137 d", "7\u201330 d", "> 30 d"))
  ) %>%
  filter(!is.na(graz_window))

nee_lim <- c(-40, 40)

p2 <- ggplot(df_wide_jc1,
             aes(x = NEE_orig_JC1, y = NEE_orig_JC2,
                 colour = Fertiliser_days_since_JC1)) +
  geom_abline(slope = 1, intercept = 0,
              colour = "grey50", linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(size = 0.85, alpha = 0.5) +
  facet_wrap(~ graz_window, nrow = 2) +
  labs(title = "Grazing at JC1") + 
  scale_colour_viridis_c(
    option    = "viridis",
    direction = -1,
    name      = "dsf at JC1"
  ) +
  scale_x_continuous(
    limits = nee_lim,
    name   = expression("JC1 NEE ("*mu*mol~m^{-2}~s^{-1}*")")
  ) +
  scale_y_continuous(
    limits = nee_lim,
    name   = expression("JC2 NEE ("*mu*mol~m^{-2}~s^{-1}*")")
  ) +
  base_theme +
  theme(legend.position   = "right",
        legend.key.height = unit(1.8, "cm"),
        legend.text       = element_text(size = 13),
        legend.title      = element_text(size = 16, face = "bold"),
        plot.title        = element_text(size = 16, face = "bold"))

save_fig(p2, "fig2_nee_jc1_vs_jc2_grazing_windows_jc1mgmt", w = 10, h = 9)


# --- FIGURE 3 — same scatter, faceted by JC2 grazing window, coloured by JC2 dsf ---

df_wide_jc2 <- df_wide_base %>%
  filter(!is.na(Grazing_days_since_JC2),
         !is.na(Fertiliser_days_since_JC2)) %>%
  mutate(
    graz_window = case_when(
      Grazing_days_since_JC2 <= 1                                  ~ "0\u20131 d",
      Grazing_days_since_JC2 >  1 & Grazing_days_since_JC2 <= 7   ~ "1\u20137 d",
      Grazing_days_since_JC2 >  7 & Grazing_days_since_JC2 <= 30  ~ "7\u201330 d",
      Grazing_days_since_JC2 >  30                                 ~ "> 30 d",
      TRUE ~ NA_character_
    ),
    graz_window = factor(graz_window,
                         levels = c("0\u20131 d", "1\u20137 d", "7\u201330 d", "> 30 d"))
  ) %>%
  filter(!is.na(graz_window))

p3 <- ggplot(df_wide_jc2,
             aes(x = NEE_orig_JC1, y = NEE_orig_JC2,
                 colour = Fertiliser_days_since_JC2)) +
  geom_abline(slope = 1, intercept = 0,
              colour = "grey50", linewidth = 0.5, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_point(size = 0.85, alpha = 0.5) +
  facet_wrap(~ graz_window, nrow = 2) +
  labs(title = "Grazing at JC2") +
  scale_colour_viridis_c(
    option    = "viridis",
    direction = -1,
    name      = "dsf at JC2"
  ) +
  scale_x_continuous(
    limits = nee_lim,
    name   = expression("JC1 NEE ("*mu*mol~m^{-2}~s^{-1}*")")
  ) +
  scale_y_continuous(
    limits = nee_lim,
    name   = expression("JC2 NEE ("*mu*mol~m^{-2}~s^{-1}*")")
  ) +
  base_theme +
  theme(legend.position   = "right",
        legend.key.height = unit(1.8, "cm"),
        legend.text       = element_text(size = 13),
        legend.title      = element_text(size = 16, face = "bold"),
        plot.title        = element_text(size = 16, face = "bold"))

save_fig(p3, "fig3_nee_jc1_vs_jc2_grazing_windows_jc2mgmt", w = 10, h = 9)


# Fig 4 — p2 and p3 side by side.
p_combined <- p2 | p3
ggsave(
  file.path(OUT_DIR, "fig4_scatter_combined.png"),
  p_combined,
  width  = 14,
  height = 6,
  dpi    = 300,
  bg     = "white",
  scale  = 1 / 0.8   # = 1.49 — elements sized so 67% zoom looks normal
)
