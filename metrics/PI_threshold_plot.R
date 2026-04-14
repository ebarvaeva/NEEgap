# =============================================================================
# PI DAYTIME PPFD THRESHOLD SENSITIVITY — PPFD > 400 vs PPFD > 700
#
# PURPOSE
#   Evaluates how sensitive the Photosynthetic Index (PI) is to the choice
#   of daytime PPFD threshold. PI is computed on the full observed NEE_orig
#   series (no gap masking) under two thresholds: PPFD > 400 and PPFD > 700
#   µmol m-2 s-1. A scatterplot of PI_400 vs PI_700 is produced for each
#   site, with a 1:1 reference line. Points above the line indicate that
#   the looser threshold (> 400) yields a higher PI estimate.
#   Both sites are shown side by side with a single shared legend.
#   Points are coloured by season (Winter, Spring, Summer, Autumn).
#
# INPUT
#   data/data_prepared/JC1_cv.rds
#   data/data_prepared/JC2_cv.rds
#
# OUTPUT
#   graphs/PI/threshold_comparison_JC1_JC2.png
# =============================================================================

library(here)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(zoo)
library(patchwork)

# ---- Helper: compute PI on whole NEE_orig series ----------------------------
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

# ---- Season colours ---------------------------------------------------------
SEASON_COLOURS <- c(
  Winter = "#4477AA",   # blue
  Spring = "#228B22",   # green
  Summer = "#FFD700",   # yellow
  Autumn = "#8B4513"    # brown
)

# ---- Helper: derive season label from dummy columns -------------------------
get_season <- function(df) {
  dplyr::case_when(
    df$Winter == 1 ~ "Winter",
    df$Spring == 1 ~ "Spring",
    df$Summer == 1 ~ "Summer",
    df$Autumn == 1 ~ "Autumn",
    TRUE           ~ NA_character_
  )
}

# ---- Helper: build scatterplot for one site (legend kept for extraction) ----
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
         title  = site_label) +                         # bold via plot.title
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "bottom",
          plot.title       = element_text(face = "bold", hjust = 0.5))
}

# ---- Build both plots, combine with patchwork -------------------------------
p_jc1 <- build_pi_plot("JC1_cv.rds", "JC1")
p_jc2 <- build_pi_plot("JC2_cv.rds", "JC2")

combined <- (p_jc1 + p_jc2 + guide_area()) +
  plot_layout(
    design  = "AB\nCC",          # A and B side by side, C spans full width below
    guides  = "collect",
    heights = c(10, 1)           # tall plots, thin legend row
  )

# ---- Save -------------------------------------------------------------------
out_dir <- here("graphs", "PI")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_path <- file.path(out_dir, "threshold_comparison_JC1_JC2.png")
ggsave(out_path, combined, width = 12, height = 6.5, dpi = 220, bg = "white")
message("Saved: ", out_path)

# =============================================================================
# end PI_threshold_comparison.R
# =============================================================================