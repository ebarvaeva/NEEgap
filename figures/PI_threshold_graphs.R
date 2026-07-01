# PI_threshold_plot.R — PI daytime-PPFD threshold sensitivity (paper Fig 3)
#
# Standalone plotting script. Recomputes the Phytomass Index (PI) on the full
# observed NEE series under two daytime PPFD thresholds (>400 and >700
# umol m-2 s-1) and draws a PI_400-vs-PI_700 scatter per site with a 1:1 line,
# so the reader can see how little the threshold choice moves PI. Both sites
# share one legend; points are coloured by season.
#
# Inputs (in data/data_prepared/):
#   JC1_cv.rds, JC2_cv.rds
#   required columns: timestamp, NEE_orig, PPFD, Winter, Spring, Summer, Autumn
#
# Output:
#   graphs/PI/threshold_comparison_JC1_JC2.png
#
# To adapt: add a site by calling build_pi_plot("<SITE>_cv.rds", "<SITE>") and
# placing it in the patchwork layout below.

library(here)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(zoo)
library(patchwork)

# compute_pi_global(): PI on the full NEE series — 21-day centred rolling mean of
# (night NEE - day NEE), scaled to [0,1]; day/night set by PPFD thresholds.
compute_pi_global <- function(df, night_ppfd = 1, day_ppfd = 400,
                              window_days = 21L) {
  
  cd  <- as.Date(df$timestamp)
  nv  <- as.numeric(df$NEE_orig)
  pv  <- as.numeric(df$PPFD)
  
  all_dates <- tibble::tibble(date = seq(min(cd), max(cd), by = "day"))
  rol <- function(x) zoo::rollapply(x, width = window_days, FUN = sum,
                                    align = "center", fill = NA_real_)
  
  ok  <- is.finite(nv) & is.finite(pv)
  isn <- ok & pv < night_ppfd
  isd <- ok & pv > day_ppfd
  
  dt <- tibble::tibble(
    date = cd,
    nn   = dplyr::if_else(isn, nv, 0), cn = as.integer(isn),
    nd   = dplyr::if_else(isd, nv, 0), cd = as.integer(isd)
  ) |>
    dplyr::group_by(date) |>
    dplyr::summarise(nn = sum(nn), cn = sum(cn),
                     nd = sum(nd), cd = sum(cd), .groups = "drop") |>
    dplyr::arrange(date)
  
  dt_full <- all_dates |>
    dplyr::left_join(dt, by = "date") |>
    dplyr::mutate(dplyr::across(c(nn, cn, nd, cd), ~ tidyr::replace_na(., 0)))
  
  mn <- dplyr::if_else(rol(dt_full$cn) > 0,
                       rol(dt_full$nn) / rol(dt_full$cn), NA_real_)
  md <- dplyr::if_else(rol(dt_full$cd) > 0,
                       rol(dt_full$nd) / rol(dt_full$cd), NA_real_)
  PR <- dplyr::if_else(is.finite(mn) & is.finite(md), mn - md, NA_real_)
  mx <- suppressWarnings(max(PR, na.rm = TRUE))
  
  PN <- if (!is.finite(mx) || mx <= 0) rep(NA_real_, length(PR)) else
    pmax(0, pmin(PR / mx, 1))
  
  PN[match(cd, dt_full$date)]
}

# Season -> colour lookup.
SEASON_COLOURS <- c(
  Winter = "#4477AA",   # blue
  Spring = "#228B22",   # green
  Summer = "#FFD700",   # yellow
  Autumn = "#8B4513"    # brown
)

# get_season(): collapse the four season dummy columns into one label.
get_season <- function(df) {
  dplyr::case_when(
    df$Winter == 1 ~ "Winter",
    df$Spring == 1 ~ "Spring",
    df$Summer == 1 ~ "Summer",
    df$Autumn == 1 ~ "Autumn",
    TRUE           ~ NA_character_
  )
}

# build_pi_plot(): load one site, compute PI at both thresholds, return the scatter.
build_pi_plot <- function(site_file, site_label) {
  
  df <- readRDS(here("data/data_prepared", site_file)) |>
    filter(!is.na(NEE_orig))
  
  df_compare <- tibble::tibble(
    PI_400 = compute_pi_global(df, day_ppfd = 400),
    PI_700 = compute_pi_global(df, day_ppfd = 700),
    season = factor(get_season(df),
                    levels = c("Winter", "Spring", "Summer", "Autumn"))
  ) |>
    filter(is.finite(PI_400) & is.finite(PI_700))
  
  ggplot(df_compare, aes(x = PI_700, y = PI_400, colour = season)) +
    geom_point(alpha = 0.3, size = 0.8) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    scale_colour_manual(values = SEASON_COLOURS, drop = FALSE,
                        guide  = guide_legend(override.aes = list(alpha = 1,
                                                                  size  = 2))) +
    labs(x      = "PI (PPFD > 700)",
         y      = "PI (PPFD > 400)",
         colour = "Season",
         title  = site_label) +
    theme_bw(base_size = 15) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "bottom",
          legend.text      = element_text(size = 14),   
          legend.title     = element_text(size = 14),   
          plot.title       = element_text(face = "bold", hjust = 0.5))
}

# One scatter per site.
p_jc1 <- build_pi_plot("JC1_cv.rds", "JC1")
p_jc2 <- build_pi_plot("JC2_cv.rds", "JC2")

# Side-by-side with a single shared legend spanning the bottom.
combined <- (p_jc1 + p_jc2 + guide_area()) +
  plot_layout(
    design  = "AB\nCC",          # A and B side by side, C spans full width below
    guides  = "collect",
    heights = c(10, 1)           # tall plots, thin legend row
  )

# Save.
out_dir <- here("graphs", "PI")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_path <- file.path(out_dir, "threshold_comparison_JC1_JC2.png")
ggsave(out_path, combined, width = 12, height = 6.5, dpi = 220, bg = "white")
