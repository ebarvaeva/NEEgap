# =============================================================================
# plot_site_comparison.R — JC1 vs JC2 Comparative Figures for Paper
# =============================================================================
#
# Produces 3 publication-ready figures comparing the two grassland sites:
#
#   Fig 1 — NEE time series with annotated management events (one panel per site)
#   Fig 2 — JC1 vs JC2 NEE scatter plots coloured by days since fertiliser,
#            faceted into 2×2 panels by grazing recovery window (0–1, 1–7,
#            7–30, > 30 d); both management variables taken from JC1 only
#   Fig 3 — Same as Fig 2 but both management variables taken from JC2 only
#
# Saved to: graphs/nee_sites/
#
# HOW TO RUN
#   Set DATA_DIR to the folder containing JC1.rds and JC2.rds, then:
#   source("scripts/plot_site_comparison.R")
# =============================================================================


# =============================================================================
# SECTION 1 — Settings
# =============================================================================

DATA_DIR <- here::here("data", "data_prepared")
OUT_DIR  <- here::here("graphs", "nee_sites")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SITES      <- c("JC1", "JC2")
SITE_COLS  <- c(JC1 = "#377EB8", JC2 = "#E41A1C")   # blue / red

# Management window: rows within N days of an event are "post-event"
POST_GRAZ_DAYS <- 30
POST_FERT_DAYS <- 14

# Light-response curve fitting range
PPFD_MAX <- 1800

# Season order (Irish meteorological)
SEASON_LEVELS <- c("Winter", "Spring", "Summer", "Autumn")


# =============================================================================
# SECTION 2 — Packages
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(lubridate); library(scales); library(patchwork)
})

save_fig <- function(p, name, w = 14, h = 7) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p,
         width = w, height = h)
  message("Saved: ", name)
}

base_theme <- theme_bw(base_size = 16) +   # raised from 12
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


# =============================================================================
# SECTION 3 — Load and prepare data
# =============================================================================

df_all <- map_dfr(SITES, function(s) {
  path <- file.path(DATA_DIR, paste0(s, ".rds"))
  if (!file.exists(path)) { message("Not found: ", path); return(NULL) }
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

message("Loaded: ", nrow(df_all), " rows | sites: ",
        paste(levels(df_all$site), collapse = ", "))


# =============================================================================
# FIGURE 1 — NEE time series with management event annotations
# =============================================================================

message("\n[Fig 1] NEE time series ...")

fig1_panels <- map(SITES, function(s) {
  d <- df_all %>% filter(site == s)
  
  ts_graz <- d %>%
    filter(!is.na(Grazing_days_since), Grazing_days_since == 0) %>%
    pull(timestamp)
  ts_fert <- d %>%
    filter(!is.na(Fertiliser_days_since), Fertiliser_days_since == 0) %>%
    pull(timestamp)
  
  ggplot(d %>% filter(!is.na(NEE_orig)),
         aes(x = timestamp, y = NEE_orig)) +
    geom_vline(xintercept = as.numeric(ts_graz),
               colour = "#4DAF4A", linewidth = 0.35, linetype = "dashed", alpha = 0.8) +
    geom_vline(xintercept = as.numeric(ts_fert),
               colour = "#A65628", linewidth = 0.85, linetype = "dashed", alpha = 1) +
    geom_point(colour = "grey20", size = 0.3, alpha = 0.4) +
    geom_vline(data = data.frame(
      xint = as.POSIXct(c(NA, NA)),
      Management = factor(c("Grazing", "Fertiliser"),
                          levels = c("Grazing", "Fertiliser"))),
      aes(xintercept = xint, colour = Management),
      linetype = "dashed", linewidth = 0.8,
      inherit.aes = FALSE, na.rm = TRUE) +
    scale_colour_manual(values = c(Grazing = "#4DAF4A", Fertiliser = "#A65628"),
                        name = "Management") +
    scale_x_datetime(
      date_breaks = "1 month",
      limits = c(
        as.POSIXct(format(min(df_all$timestamp, na.rm = TRUE), "%Y-01-01"),
                   tz = "UTC"),
        max(df_all$timestamp, na.rm = TRUE)
      ),
      expand = expansion(mult = 0.01),
      labels = function(x) {
        ifelse(format(x, "%m") == "01",
               paste0(format(x, "%b"), "\n", format(x, "%Y")),
               format(x, "%b"))
      }
    ) +
    scale_y_continuous(limits = c(-40, 40),
                       name = expression(NEE~(mu*mol~m^{-2}~s^{-1}))) +
    labs(title = s, x = NULL) +
    base_theme +
    theme(legend.position = if (s == SITES[length(SITES)]) "right" else "none",
          plot.title = element_text(face = "bold", hjust = 0))
})

p1 <- wrap_plots(fig1_panels, ncol = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(theme = theme(plot.title = element_blank()))

save_fig(p1, "fig1_nee_timeseries", w = 15, h = 9)


# =============================================================================
# FIGURE 2 — JC1 vs JC2 NEE scatter coloured by days since fertiliser,
#            faceted into 2×2 panels by grazing recovery window.
#
#   Panel facets : Grazing_days_since_JC1    (rows where JC1 is NA are dropped)
#   Point colour : Fertiliser_days_since_JC1 (rows where JC1 is NA are dropped)
# =============================================================================

message("\n[Fig 2] JC1 vs JC2 NEE scatter by grazing window (JC1 management) ...")

df_wide_base <- df_all %>%
  filter(!is.na(NEE_orig)) %>%
  select(timestamp, site, NEE_orig, Grazing_days_since, Fertiliser_days_since) %>%
  pivot_wider(
    id_cols     = timestamp,
    names_from  = site,
    values_from = c(NEE_orig, Grazing_days_since, Fertiliser_days_since)
  ) %>%
  filter(!is.na(NEE_orig_JC1), !is.na(NEE_orig_JC2))

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


# =============================================================================
# FIGURE 3 — JC1 vs JC2 NEE scatter coloured by days since fertiliser,
#            faceted into 2×2 panels by grazing recovery window.
#
#   Panel facets : Grazing_days_since_JC2    (rows where JC2 is NA are dropped)
#   Point colour : Fertiliser_days_since_JC2 (rows where JC2 is NA are dropped)
# =============================================================================

message("\n[Fig 3] JC1 vs JC2 NEE scatter by grazing window (JC2 management) ...")

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


# Optional: p2 and p3 side-by-side
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
# =============================================================================
# Summary
# =============================================================================

message("\n", strrep("=", 60))
message("All figures saved to: ", OUT_DIR)
message(strrep("=", 60))
message("  fig1_nee_timeseries.{png,pdf}                          — NEE time series + management annotations")
message("  fig2_nee_jc1_vs_jc2_grazing_windows_jc1mgmt.{png,pdf} — JC1 vs JC2 NEE scatter, JC1 management")
message("  fig3_nee_jc1_vs_jc2_grazing_windows_jc2mgmt.{png,pdf} — JC1 vs JC2 NEE scatter, JC2 management")

# ==================== end plot_site_comparison.R ==============================