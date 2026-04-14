# =============================================================================
# plot_real_gap_cleveland.R — Real Gap Counts: Cleveland Dot Plot JC1 vs JC2
# =============================================================================
#
# Y-axis: gap length in 1-day bins from 0–1 d up to 29–30 d, plus a ">30 d"
#         catch-all — giving a continuous 0–30 d view of the distribution.
# X-axis: number of gaps (log scale).
#
# Output: graphs/real_gaps/histogram/real_gap_cleveland.{png,pdf}
# =============================================================================


# =============================================================================
# SECTION 1 — Settings
# =============================================================================

DATA_DIR   <- here::here("data", "data_prepared")
OUT_DIR    <- here::here("graphs", "real_gaps", "histogram")

SITES      <- c("JC1", "JC2")
SITE_COLS  <- c(JC1 = "#377EB8", JC2 = "#E41A1C")
SITE_SHAPE <- c(JC1 = 21, JC2 = 23)   # circle / diamond

HH_PER_DAY <- 48   # half-hours per day
MAX_DAYS   <- 30   # bins run 0–1, 1–2, …, 29–30, then >30


# =============================================================================
# SECTION 2 — Packages & helpers
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(lubridate)
})

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(p, name, w = 8, h = 10) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p,
         width = w, height = h)
  message("Saved: ", name)
}

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88"),
    plot.background    = element_rect(fill = "white", colour = NA),
    legend.position    = "bottom",
    axis.ticks.y       = element_blank()
  )


# =============================================================================
# SECTION 3 — Build 1-day bins
# =============================================================================

# Breaks in half-hours: 0, 48, 96, …, 1440, Inf
bin_breaks <- c(seq(0, MAX_DAYS * HH_PER_DAY, by = HH_PER_DAY), Inf)

# Labels: "0–1 d", "1–2 d", …, "29–30 d", ">30 d"
bin_labels <- c(
  paste0(0:(MAX_DAYS - 1), "\u2013", 1:MAX_DAYS, " d"),
  ">30 d"
)

gap_df <- map_dfr(SITES, function(s) {
  path <- file.path(DATA_DIR, paste0(s, ".rds"))
  if (!file.exists(path)) { message("Not found: ", path); return(NULL) }
  df <- readRDS(path) %>% arrange(as.POSIXct(timestamp, tz = "UTC"))
  r  <- rle(is.na(df$NEE_orig))
  tibble(site = s, length_hh = r$lengths[r$values]) %>%
    mutate(category = cut(length_hh,
                          breaks         = bin_breaks,
                          labels         = bin_labels,
                          include.lowest = TRUE,
                          right          = TRUE))
}) %>%
  mutate(site = factor(site, levels = SITES))

# Complete all combinations — zeros explicit; reverse so 0-1d is at top
plot_df <- gap_df %>%
  count(site, category) %>%
  complete(site,
           category = factor(bin_labels, levels = bin_labels),
           fill     = list(n = 0)) %>%
  mutate(category = factor(category, levels = rev(bin_labels))) %>%
  group_by(category) %>%
  filter(sum(n) > 0) %>%
  ungroup()

print(plot_df, n = Inf)


# =============================================================================
# SECTION 4 — Cleveland dot plot
# =============================================================================

p <- ggplot(plot_df, aes(x = n, y = category,
                         colour = site, shape = site)) +
  # connecting segment
  geom_line(aes(group = category), colour = "grey75", linewidth = 0.5) +
  # dots
  geom_point(size = 2.5, stroke = 0.9, fill = "white") +
  # count labels: only show where n > 0 to avoid clutter
  geom_text(
    data = plot_df %>% filter(n > 0, site == "JC1"),
    aes(label = n), nudge_y =  0.3, size = 2.8, show.legend = FALSE
  ) +
  geom_text(
    data = plot_df %>% filter(n > 0, site == "JC2"),
    aes(label = n), nudge_y = -0.3, size = 2.8, show.legend = FALSE
  ) +
  scale_x_log10(
    name   = "Number of gaps  (log scale)",
    labels = scales::comma,
    breaks = c(1, 10, 100, 1000, 10000),
    limits = c(0.8, 15000)
  ) +
  scale_colour_manual(values = SITE_COLS,  name = "Site") +
  scale_shape_manual( values = SITE_SHAPE, name = "Site") +
  labs(title = "Distribution of real gap lengths", y = "Gap length") +
  base_theme +
  guides(
    colour = guide_legend(override.aes = list(size = 3, stroke = 1,
                                              fill = "white")),
    shape  = guide_legend(override.aes = list(size = 3, stroke = 1,
                                              fill = "white"))
  )

save_fig(p, "real_gap_cleveland")
message("\nDone. Output: ", OUT_DIR)
# =================== end plot_real_gap_cleveland.R ===========================