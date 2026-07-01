# =============================================================================
# plot_artificial_gaps.R — Artificial Gap Positions on NEE Time Series
# =============================================================================
#
# Saves:
#   graphs/artificial_gaps/{SITE}/gap_positions_{SITE}_{SIZE}.{png,pdf}
#     — one file per site × gap-size (S, M, L, VL)
#
#   graphs/artificial_gaps/comparison/gap_positions_comparison_{SIZE}.{png,pdf}
#     — one file per gap-size: JC1 panel on top, JC2 panel below
#
# Input:  results/{SITE}/managed/RF/df_cv_all_predictions.rds
#           columns: timestamp, NEE_orig, S1..Sn, M1..Mn, L1..Ln, VL1..VLn
#           (binary 0/1: 1 = observation belongs to that gap)
#
# HOW TO RUN
#   source("scripts/plot_artificial_gaps.R")
# =============================================================================


# =============================================================================
# SECTION 1 — Settings
# =============================================================================

SITES     <- c("JC1", "JC2")
GAP_SIZES <- c("S", "M", "L", "VL")

GAP_LABELS <- c(S  = "Small (1-day)",
                M  = "Medium (7-day)",
                L  = "Large (14-day)",
                VL = "Very Large (30-day)")

NEE_LIM <- c(-40, 40)


# =============================================================================
# SECTION 2 — Packages & helpers
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(purrr)
  library(ggplot2); library(lubridate); library(patchwork)
})

save_fig <- function(p, name, dir, w = 18, h = 7) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, paste0(name, ".png")), p,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(dir, paste0(name, ".pdf")), p,
         width = w, height = h)
  message("  Saved: ", file.path(dir, name))
}

base_theme <- theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    plot.background  = element_rect(fill = "white", colour = NA),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold", size = 13),
    legend.position  = "none"
  )

x_date_scale <- function(data) {
  scale_x_datetime(
    date_breaks = "1 month",
    limits      = c(
      as.POSIXct(format(min(data$timestamp, na.rm = TRUE), "%Y-01-01"), tz = "UTC"),
      max(data$timestamp, na.rm = TRUE)
    ),
    expand  = expansion(mult = 0.01),
    labels  = function(x) {
      ifelse(format(x, "%m") == "01",
             paste0(format(x, "%b"), "\n", format(x, "%Y")),
             format(x, "%b"))
    }
  )
}


# =============================================================================
# SECTION 3 — Gap-start extractor
# =============================================================================

# For prefix "S": finds S1, S2, ..., returns first timestamp == 1 per column
gap_starts <- function(df, prefix) {
  cols <- grep(paste0("^", prefix, "[0-9]+$"), names(df), value = TRUE)
  cols <- cols[order(as.integer(sub(prefix, "", cols)))]
  map_dfr(cols, function(col) {
    idx <- which(df[[col]] == 1)
    if (length(idx) == 0) return(NULL)
    tibble(gap_id  = col,
           t_start = df$timestamp[min(idx)])
  })
}


# =============================================================================
# SECTION 4 — Single-panel builder  (one site x one size)
# =============================================================================

build_panel <- function(df, site, size, show_x = TRUE) {
  starts <- gap_starts(df, size)
  n_gaps <- nrow(starts)
  message("    ", site, " / ", size, ": ", n_gaps, " gaps")
  
  p <- ggplot(df %>% filter(!is.na(NEE_orig)),
              aes(x = timestamp, y = NEE_orig)) +
    geom_vline(data        = starts,
               aes(xintercept = as.numeric(t_start)),
               colour      = "black",
               linewidth   = 0.35,
               alpha       = 0.55,
               inherit.aes = FALSE) +
    geom_point(colour = "grey30", size = 0.25, alpha = 0.35) +
    x_date_scale(df) +
    scale_y_continuous(limits = NEE_LIM,
                       name   = expression(NEE~(mu*mol~m^{-2}~s^{-1}))) +
    labs(title = paste0(site, "  —  ", GAP_LABELS[size],
                        "  (n\u202f=\u202f", n_gaps, " gaps)"),
         x     = NULL) +
    base_theme +
    theme(plot.title = element_text(face = "bold", hjust = 0, size = 15))
  
  # suppress x-axis text/ticks on upper panels in stacked comparisons
  if (!show_x)
    p <- p + theme(axis.text.x  = element_blank(),
                   axis.ticks.x = element_blank())
  p
}


# =============================================================================
# SECTION 5 — Load data for both sites
# =============================================================================

df_list <- map(SITES, function(s) {
  path <- here::here("results", s, "managed", "RF",
                     "df_cv_all_predictions.rds")
  if (!file.exists(path)) {
    message("Not found: ", path, " — skipping ", s)
    return(NULL)
  }
  readRDS(path) %>%
    mutate(timestamp = as.POSIXct(timestamp, tz = "UTC"))
}) %>% set_names(SITES)


# =============================================================================
# SECTION 6 — Individual plots  (one file per site x size)
# =============================================================================

message("\n── Individual plots ──")

walk(SITES, function(s) {
  if (is.null(df_list[[s]])) return(invisible())
  out_dir <- here::here("graphs", "artificial_gaps", s)
  walk(GAP_SIZES, function(size) {
    p <- build_panel(df_list[[s]], s, size)
    save_fig(p,
             paste0("gap_positions_", s, "_", size),
             out_dir,
             w = 18, h = 5)
  })
})


# =============================================================================
# SECTION 7 — Comparison plots  (JC1 top / JC2 bottom, one file per size)
# =============================================================================

message("\n── Comparison plots ──")

comp_dir <- here::here("graphs", "artificial_gaps", "comparison")

walk(GAP_SIZES, function(size) {
  panels <- map(SITES, function(s) {
    if (is.null(df_list[[s]])) return(NULL)
    # only bottom panel (last site) gets x-axis labels
    show_x <- (s == SITES[length(SITES)])
    build_panel(df_list[[s]], s, size, show_x = show_x)
  }) %>% compact()
  
  if (length(panels) == 0) return(invisible())
  
  p_comp <- wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title = paste0("Artificial gap positions — ", GAP_LABELS[size]),
      theme = theme(plot.title = element_text(face = "bold",
                                              size = 15, hjust = 0))
    )
  
  save_fig(p_comp,
           paste0("gap_positions_comparison_", size),
           comp_dir,
           w = 18, h = 10)
})

message("\nDone.")
# ========================= end plot_artificial_gaps.R =========================