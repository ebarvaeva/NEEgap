# wind_rose_graphs.R — Wind roses per site and year (JC1 2020/2023/2024, JC2 2023/2024)
#
# Loads QC'd half-hourly wind speed and direction for each site x year, bins the
# data into 16 direction sectors x 4 speed classes, and draws a faceted wind
# rose (site rows, year columns). Rings, compass letters, ring labels and the
# per-panel mean are all drawn as data so their radii can be tuned; each panel's
# wedges sum to 100%.
#
# Input:  data/data_qc/results_qc_JC{site}_{year}/merged_qc.rds
#           (needs wind_speed, wind_dir, and date or timestamp)
# Output: graphs/wind_rose/windrose_site_year.png
#
# To adapt: edit the load_and_prepare() calls for other site/year combinations;
# sector / brks / labs set the binning; the r_* multipliers place the overlays.

# Packages.
library(dplyr)
library(ggplot2)
library(ggh4x)      # install.packages("ggh4x")  -- render_empty = FALSE
library(scales)     # oob_keep

# Output directory.
dir.create("graphs/wind_rose", recursive = TRUE, showWarnings = FALSE)

# load_and_prepare(): read one site x year QC file and return clean wind rows.
load_and_prepare <- function(site, year) {
  
  path <- sprintf("data/data_qc/results_qc_JC%d_%d/merged_qc.rds", site, year)
  df <- readRDS(path)
  
  # Ensure a date column exists (fall back to timestamp).
  if (!"date" %in% names(df)) {
    if ("timestamp" %in% names(df)) {
      df$date <- df$timestamp
    } else {
      stop("No 'date' or 'timestamp' column found")
    }
  }
  
  # Coerce wind fields, wrap direction to 0-359, drop bad/missing rows.
  df %>%
    mutate(
      wind_speed = as.numeric(wind_speed),
      wind_dir   = as.numeric(wind_dir) %% 360,   # 360 -> 0, no double-counted N
      site       = paste0("JC", site),
      year       = as.character(year)
    ) %>%
    filter(
      !is.na(wind_speed),
      !is.na(wind_dir),
      wind_speed >= 0
    )
}

# Load every site x year and lock factor order.
data_all <- bind_rows(
  load_and_prepare(1, 2020),
  load_and_prepare(1, 2023),
  load_and_prepare(1, 2024),
  load_and_prepare(2, 2023),
  load_and_prepare(2, 2024)
) %>%
  mutate(
    site = factor(site, levels = c("JC1", "JC2")),
    year = factor(year, levels = c("2020", "2023", "2024"))
  )

# Bin into 16 direction sectors (22.5 deg each) x 4 speed classes.
sector <- 22.5
brks   <- c(0, 2, 4, 6, Inf)
labs   <- c("0-2", "2-4", "4-6", ">6")

# Count sector x speed bins, then convert to per-panel % frequency.
rose <- data_all %>%
  mutate(
    wd_bin = (round(wind_dir / sector) * sector) %% 360,
    ws_bin = cut(wind_speed, breaks = brks, labels = labs,
                 right = FALSE, include.lowest = TRUE)
  ) %>%
  count(site, year, wd_bin, ws_bin, name = "n") %>%
  group_by(site, year) %>%
  mutate(freq = 100 * n / sum(n)) %>%          # each panel sums to 100%
  ungroup()

# Radial extent, rings, and label radii.
# Multipliers act on the radial span: 1.00 == the outermost ring.
rmax <- rose %>%
  group_by(site, year, wd_bin) %>%
  summarise(tot = sum(freq), .groups = "drop") %>%
  pull(tot) %>% max()
rmax <- ceiling(rmax / 5) * 5

rmin  <- -2                            # the centre hole
rings <- seq(0, rmax, by = 5)

r_dir <- rmin + 1.18 * (rmax - rmin)   # compass letters (1.12 ~ ggplot default)
r_lab <- rmin + 1.55 * (rmax - rmin)   # mean, lower-right corner

# Ring labels, drawn along the NW ray in every panel.
ring_lab <- data.frame(y = rings, label = paste0(rings, "%"))

# Compass letters, drawn as data so their radius is controllable.
compass <- data.frame(
  x     = seq(0, 315, 45),
  label = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")
)

# Per-panel mean; sprintf, not round(): round(2.9046, 2) prints "2.9", not "2.90".
mean_lab <- data_all %>%
  group_by(site, year) %>%
  summarise(mean_ws = mean(wind_speed), .groups = "drop") %>%
  mutate(label = sprintf("mean = %.2f", mean_ws))

# Full summary for the caption.
summary_tbl <- data_all %>%
  group_by(site, year) %>%
  summarise(
    n        = n(),
    mean_ws  = sprintf("%.2f", mean(wind_speed)),
    calm_pct = sprintf("%.2f", 100 * mean(wind_speed < 0.5)),
    .groups  = "drop"
  )
print(summary_tbl)

# Build the faceted wind rose.
p <- ggplot(rose, aes(x = wd_bin, y = freq, fill = ws_bin)) +
  
  # rings drawn as data, first, so wedges sit on top of them
  geom_hline(yintercept = rings, colour = "grey88", linewidth = 0.3) +
  
  geom_col(
    width     = sector,
    colour    = "white",
    linewidth = 0.2,
    position  = position_stack(reverse = TRUE)
  ) +
  
  # ring labels along the NW ray
  geom_label(
    data          = ring_lab,
    mapping       = aes(x = 315, y = y, label = label),
    inherit.aes   = FALSE,
    size          = 2.3,
    colour        = "grey40",
    fill          = "white",
    alpha         = 0.7,
    label.size    = 0,
    label.padding = unit(0.06, "lines")
  ) +
  
  # compass letters, pushed out past the outer ring
  geom_text(
    data        = compass,
    mapping     = aes(x = x, y = r_dir, label = label),
    inherit.aes = FALSE,
    size        = 3.5,                 # mm, not pt: 4.4 mm ~ 12.4 pt
    colour      = "grey20"
  ) +
  
  # mean wind speed, lower right, nudged off the SE ray
  geom_text(
    data        = mean_lab,
    mapping     = aes(x = 157.5, y = r_lab, label = label),
    inherit.aes = FALSE,
    size        = 2.8,
    colour      = "grey20"
  ) +
  
  # x = 0 (N) rotated back to 12 o'clock, since lower limit is -11.25
  coord_polar(theta = "x", start = -sector * pi / 180 / 2, clip = "off") +
  
  scale_x_continuous(
    limits = c(-sector / 2, 360 - sector / 2),
    breaks = seq(0, 315, 45)
  ) +
  scale_y_continuous(
    limits = c(rmin, rmax),
    breaks = rings,
    expand = c(0, 0),                  # no phantom ring at the panel edge
    oob    = scales::oob_keep          # without this, labels past rmax vanish
  ) +
  
  scale_fill_viridis_d(
    option    = "viridis",
    direction = -1,                    # dark at the rim, where wedges meet white
    name      = expression(paste("Wind speed (m s"^-1, ")"))
  ) +
  
  # blank JC2 x 2020 cell: not rendered at all, no empty circle
  ggh4x::facet_grid2(site ~ year, render_empty = FALSE) +
  
  guides(fill = guide_legend(
    nrow           = 1,
    title.position = "left",
    label.position = "bottom",
    keywidth       = unit(2.6, "lines")
  )) +
  
  theme_minimal(base_size = 13) +
  theme(
    axis.title         = element_blank(),
    axis.text.x        = element_blank(),   # compass drawn as data instead
    axis.text.y        = element_blank(),   # ring labels drawn as data instead
    panel.grid.major.y = element_blank(),   # rings drawn as data instead
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    strip.text         = element_text(size = rel(1.05)),
    
    legend.position    = "bottom",
    legend.direction   = "horizontal",
    legend.title       = element_text(size = 10, vjust = 0.85),
    legend.key.height  = unit(0.6, "lines"),
    legend.box.spacing = unit(30, "pt"),
    legend.margin      = margin(0, 0, 0, 0),
    
    panel.spacing = unit(1.0, "lines"),     # room for compass + mean overflow
    plot.margin   = margin(6, 6, 6, 6)
  )

# Save at the printed size -- do NOT downscale in LaTeX.
ggsave("graphs/wind_rose/windrose_site_year.png", p,
       width = 8, height = 5.8, dpi = 300, bg = "white")
