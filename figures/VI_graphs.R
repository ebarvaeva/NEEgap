# VI_graphs.R — RF variable-importance summary boxplots + CSV summaries
#
# Reads pre-computed RF variable-importance RDS files (no model re-fitting) and, for
# every retained predictor, draws the full distribution of relative importance pooled
# across all artificial gaps and gap sizes (S, M, L, VL), one boxplot per site
# (JC1 blue, JC2 red). Two figure families: Plot 1 for the managed predictor set, and
# Plot 2 for each BASE + one-management-variable feature-addition run. Each figure is
# saved as PNG + PDF; matching median/Q25/Q75 summaries are also written as CSVs.
#
# Input RDS (one per gap size, columns predictor + importance_rel summing to 1 per gap):
#   Plot 1:  results/{SITE}/managed/RF/variable_importance/rf_variable_importance_{S,M,L,VL}_all.rds
#   Plot 2:  results/{SITE}/management_effect/RF/{MGMT_VAR}/variable_importance/rf_variable_importance_{S,M,L,VL}_all.rds
#   Per-gap PI columns (PI_S1, PI_M3, …) are collapsed to the single label "PI".
#
# Predictors kept if median importance_rel > IMP_THRESHOLD OR the predictor is one of
# the 5 management variables (dsg, dsf, h[s], B, PI — always shown). y-axis symbols
# come from PREDICTOR_LABELS (R plotmath, so subscripts render), ordered by descending
# JC2 median importance.
#
# Outputs:
#   graphs/VI/plot1_managed_RF.{png,pdf}
#   graphs/VI/plot2_{MGMT_VAR}_RF.{png,pdf}   — one per management variable
#   graphs/VI/csv/vi_summary_managed.csv
#   graphs/VI/csv/vi_summary_{MGMT_VAR}.csv   — one per management variable
#
# To adapt: change SITES / GAP_SIZES / MGMT_VARS (and the label lookups), or
# IMP_THRESHOLD. Run with: source("scripts/VI_graphs_compact.R")


# Output directory for the figures
RESULTS_ROOT <- here::here("results")
OUT_DIR      <- here::here("graphs", "VI")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SITES     <- c("JC1", "JC2")
GAP_SIZES <- c("S", "M", "L", "VL")

# Predictors with median importance_rel <= this are dropped unless they are one
# of the 5 management variables (which are always shown)
IMP_THRESHOLD <- 0

# Management variables for Plot 2 — order controls plot production sequence
MGMT_VARS <- c(
  "Grazing_days_since",
  "Fertiliser_days_since",
  "grass_height",
  "grass_biomass",
  "phytomass_index"
)

# Always-show list, using the unified "PI" label produced by tidy_vi_raw()
MGMT_PLOT_LABELS <- c(
  "Grazing_days_since",
  "Fertiliser_days_since",
  "grass_height",
  "grass_biomass",
  "PI"
)

# Short display labels used in Plot 2 subtitles
MGMT_COL_LABELS <- c(
  Grazing_days_since    = "Grazing",
  Fertiliser_days_since = "Fertiliser",
  grass_height          = "Height",
  grass_biomass         = "Biomass",
  phytomass_index       = "PI"
)

# Code name -> display symbol (R plotmath; X[y] = X subscript y); unlisted names shown as-is
PREDICTOR_LABELS <- c(
  # Meteorological and radiation drivers
  PPFD              = "PPFD",
  night             = "night",
  ustar             = "u['*']",          # u*
  Temp              = "T[air]",          # T_air
  RH                = "RH",
  rain              = "P",               # P
  rain_rolling_24   = "P[24]",           # P_24
  Rg                = "R[g]",            # R_g
  VPD               = "VPD",
  # Temporal cyclic encodings
  hour_sin          = "hour_sin",
  hour_cos          = "hour_cos",
  doy_sin           = "doy_sin",
  doy_cos           = "doy_cos",
  month_sin         = "month_sin",
  month_cos         = "month_cos",
  # Seasonal indicators
  Winter            = "Winter",
  Spring            = "Spring",
  Summer            = "Summer",
  Autumn            = "Autumn",
  # Management variables
  Grazing_days_since    = "bold(dsg)",
  Fertiliser_days_since = "bold(dsf)",
  grass_height          = "bold(h[s])",
  grass_biomass         = "bold(B)",
  PI                    = "bold(PI)"
)


# Packages
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(forcats); library(stringr); library(glue)
  library(scales)
})

# Save one figure as both PNG (300 dpi) and PDF
save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p,
         width = w, height = h)
}

# Read the four gap-size RDS files from a folder, returning every per-gap row (not averaged)
load_rf_vi_raw <- function(vi_root) {
  map_dfr(GAP_SIZES, function(gs) {
    path <- file.path(vi_root, "variable_importance",
                      glue("rf_variable_importance_{gs}_all.rds"))
    if (!file.exists(path)) {
      return(NULL)
    }
    readRDS(path) %>%
      mutate(gap_size = gs)
  })
}

# Collapse all PI-related predictor names (PI_*, phytomass_index) to the single label "PI"
tidy_vi_raw <- function(df) {
  df %>%
    mutate(predictor = case_when(
      str_starts(predictor, "PI_")   ~ "PI",
      predictor == "phytomass_index" ~ "PI",
      TRUE                           ~ predictor
    )) %>%
    mutate(gap_size = factor(gap_size, levels = GAP_SIZES))
}

# Keep a predictor if its median importance exceeds IMP_THRESHOLD or it is a management variable
filter_predictors <- function(df) {
  med_imp <- df %>%
    group_by(predictor) %>%
    summarise(med = median(importance_rel, na.rm = TRUE), .groups = "drop")

  keep <- med_imp %>%
    filter(med > IMP_THRESHOLD | predictor %in% MGMT_PLOT_LABELS) %>%
    pull(predictor)

  df %>% filter(predictor %in% keep)
}

# Shared boxplot builder for Plot 1 and Plot 2; x labels parsed from PREDICTOR_LABELS
make_vi_boxplot <- function(df, pred_order, subtitle = NULL) {

  df <- df %>%
    mutate(
      predictor = factor(predictor, levels = pred_order),
      site      = factor(site,      levels = SITES),
      gap_size  = factor(gap_size,  levels = GAP_SIZES)
    )

  # Map code names to plotmath symbols, falling back to the raw name when unlisted
  x_labels <- function(x) {
    labs <- PREDICTOR_LABELS[x]
    labs[is.na(labs)] <- x[is.na(labs)]
    parse(text = labs)
  }

  SITE_COLOURS <- c(JC1 = "blue", JC2 = "red")

  p <- ggplot(df,
              aes(x    = predictor,
                  y    = importance_rel,
                  fill = site)) +             # ← fill encodes site
    geom_boxplot(
      outlier.size  = 0.9,
      outlier.alpha = 0.9,
      linewidth     = 0.3,
      position      = position_dodge(width = 0.9),   # ← side-by-side
      width         = 0.9
    ) +
    coord_flip() +
    scale_fill_manual(values = SITE_COLOURS, name = "Site") +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(
      name   = "Relative importance",
      breaks = pretty_breaks(n = 6),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    labs(x = NULL, subtitle = subtitle) +
    theme_bw(base_size = 35) +
    theme(
      legend.position    = "bottom",
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background   = element_rect(fill = "grey92"),
      strip.text         = element_text(size = 32, face = "bold"),
      axis.text.x        = element_text(angle = 32, hjust = 1, size = 32),
      axis.text.y        = element_text(size = 32),
      axis.title.y       = element_text(size = 32),
      legend.text        = element_text(size = 32),
      legend.title       = element_text(size = 32),
      plot.subtitle      = element_text(size = 32, face = "bold")
    )
  p
}


# Plot 1 — managed predictor set: load both sites, tidy, and stack
vi_managed_raw <- map_dfr(SITES, function(s) {
  vi_root <- file.path(RESULTS_ROOT, s, "managed", "RF")
  raw     <- load_rf_vi_raw(vi_root)
  if (is.null(raw) || nrow(raw) == 0) {
    return(NULL)
  }
  tidy_vi_raw(raw) %>% mutate(site = s)
})

if (nrow(vi_managed_raw) == 0) stop("No managed VI data found.")

vi_managed_filtered <- filter_predictors(vi_managed_raw)

# Order predictors by ascending JC2 median (coord_flip puts the largest at the top)
predictor_order_managed <- vi_managed_filtered %>%
  filter(site == "JC2") %>%                                    # ← JC2 only
  group_by(predictor) %>%
  summarise(med = median(importance_rel, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(predictor)

p1 <- make_vi_boxplot(vi_managed_filtered, predictor_order_managed,
                      subtitle = "")

save_fig(p1, "plot1_managed_RF", w = 30, h = 17)


# Plot 2 — one figure per management-effect ablation run
for (mv in MGMT_VARS) {

  vi_mv_raw <- map_dfr(SITES, function(s) {
    vi_root <- file.path(RESULTS_ROOT, s, "management_effect", "RF", mv)
    raw     <- load_rf_vi_raw(vi_root)
    if (is.null(raw) || nrow(raw) == 0) {
      return(NULL)
    }
    tidy_vi_raw(raw) %>% mutate(site = s)
  })

  if (nrow(vi_mv_raw) == 0) {
    next
  }

  vi_mv_filtered <- filter_predictors(vi_mv_raw)

  # Order predictors by ascending JC2 median within this run
  pred_order_mv <- vi_mv_filtered %>%
    filter(site == "JC2") %>%                                    # ← JC2 only
    group_by(predictor) %>%
    summarise(med = median(importance_rel, na.rm = TRUE), .groups = "drop") %>%
    arrange(med) %>%
    pull(predictor)

  p2 <- make_vi_boxplot(
    df         = vi_mv_filtered,
    pred_order = pred_order_mv,
    subtitle   = paste0("BASE + ", MGMT_COL_LABELS[mv])
  )

  save_fig(p2, name = paste0("plot2_", mv, "_RF"),  w = 30, h = 17)

}


# CSV summaries (median, Q25, Q75 per predictor × site)
OUT_CSV <- here::here("graphs", "VI", "csv")
dir.create(OUT_CSV, recursive = TRUE, showWarnings = FALSE)

# Per predictor × site median and quartiles, sorted by descending median
compute_vi_summary <- function(df) {
  df %>%
    group_by(predictor, site) %>%
    summarise(
      median = median(importance_rel, na.rm = TRUE),
      Q25    = quantile(importance_rel, 0.25, na.rm = TRUE),
      Q75    = quantile(importance_rel, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(site, desc(median))
}

# Managed predictor set summary
compute_vi_summary(vi_managed_filtered) %>%
  write.csv(file.path(OUT_CSV, "vi_summary_managed.csv"), row.names = FALSE)

# One summary CSV per feature-addition run
for (mv in MGMT_VARS) {

  vi_mv_raw <- map_dfr(SITES, function(s) {
    vi_root <- file.path(RESULTS_ROOT, s, "management_effect", "RF", mv)
    raw     <- load_rf_vi_raw(vi_root)
    if (is.null(raw) || nrow(raw) == 0) return(NULL)
    tidy_vi_raw(raw) %>% mutate(site = s)
  })

  if (nrow(vi_mv_raw) == 0) {
    next
  }

  fname <- paste0("vi_summary_", mv, ".csv")
  filter_predictors(vi_mv_raw) %>%
    compute_vi_summary() %>%
    write.csv(file.path(OUT_CSV, fname), row.names = FALSE)
}
