# =============================================================================
# management_effect_graphs.R — Management Variable Contribution Analysis
# =============================================================================
#
# PURPOSE
#   Quantifies and visualises how much each management variable improves
#   gap-filling performance relative to the BASE (meteorological only)
#   predictor set, across three models, two sites, four gap sizes, and
#   three flux targets.
#
# STUDY DESIGN
#   BASE predictors:  PPFD, Rg, VPD, RH, Temp, rain, rain_rolling_24,
#                     temporal encodings, season dummies, night flag.
#   For each management variable, the model is retrained with BASE + one
#   additional variable.  The improvement relative to BASE is the management
#   variable's marginal contribution.  Phytomass Index (PI) is included as
#   a derived management variable alongside the five raw ones.
#
# DELTA CONVENTION
#   delta = metric(BASE + mgmt) - metric(BASE)
#   For MAE and RMSE: negative delta = improvement (lower error).
#   For R²:          positive delta = improvement (higher fit).
#   pct_mae_reduction = -100 * delta_mae / base_mae
#     (positive = improvement, printed in heatmaps).
#
# GRAZING WINDOW SPLIT
#   All metrics and deltas are computed separately for two temporal windows:
#     le30 — gap rows where Grazing_days_since <= SPLIT_DAYS (early recovery)
#     gt30 — gap rows where Grazing_days_since >  SPLIT_DAYS (late recovery)
#   Management variables are expected to show stronger improvements in the
#   early recovery window where the canopy is actively changing.
#
# INPUTS
#   graphs/metrics_csv/gap_metrics_{TARGET}_{GAP_SIZE}.csv
#     BASE (unmanaged) metrics — produced by compute_metrics_and_plots.R.
#     Run that script first; this one will stop if the CSVs are absent.
#
#   results/{SITE}/management_effect/{MODEL}/{MGMT_VAR}/df_cv_all_predictions.rds
#     Ablation run predictions — produced by the run_management_effect_*.R scripts.
#
# OUTPUTS  (saved under graphs/management_effect/)
#   delta_barplots/        — ΔMAE/ΔRMSE/ΔR² per variable, grouped bars by model
#   connected_dots/        — absolute metric; lines = predictor set; x = model
#   heatmaps/              — % MAE reduction matrix; rows = variable, cols = model
#   recovery_interaction/  — delta vs grazing recovery phase (High/Mid/Low)
#   seasonal_breakdown/    — ΔMAE boxplots by dominant season of the gap
#   ranking_stability/     — bump charts: management variable rank across models
#
# SECTIONS
#   1  Settings         — sites, models, targets, gap sizes, management vars
#   2  Libraries        — packages, metric helpers, theme, colour palettes
#   3  Load BASE        — read unmanaged metrics from CSVs
#   4  Compute metrics  — read ablation RDS files; compute le30/gt30 metrics
#   5  Build BASE join  — reformat BASE metrics for delta computation
#   6  Aggregate/delta  — weighted means per group; compute all delta columns
#   7  Plot 1           — delta barplots (ΔMAE, ΔRMSE, ΔR²)
#   8  Plot 2           — connected-dot plots (absolute metric)
#   9  Plot 3           — heatmaps (% MAE reduction)
#  10  Plot 4           — recovery-phase interaction lines
#  11  Plot 5           — seasonal breakdown boxplots
#  12  Plot 6           — ranking stability bump charts
#  13  Summary message
#
# HOW TO RUN
#   source("metrics/management_effect_graphs.R")
#   compute_metrics_and_plots.R must have been run first.
#
# =============================================================================


# =============================================================================
# SECTION 1 — Settings
# =============================================================================
# Edit the variables in this block to adapt the script.

RESULTS_ROOT <- here::here("results")
GRAPHS_ROOT  <- here::here("graphs")
OUT_ROOT     <- file.path(GRAPHS_ROOT, "management_effect")

SITES       <- c("JC1", "JC2")
MODELS      <- c("RF", "MLP", "XGBoost")
TARGET_VARS <- c("NEE", "Reco", "GPP")
GAP_SIZES   <- c("S", "M", "L", "VL")
SPLIT_DAYS  <- 30    # days since grazing defining the le30/gt30 window boundary
GRAZ_HIGH   <- 0.80  # fraction threshold for "High" recovery group
GRAZ_MID    <- 0.50  # fraction threshold for "Mid"  recovery group

# PI (phytomass_index) added as the last entry — stored under its own subfolder
MGMT_VARS <- c(
  "Grazing_days_since",
  "Fertiliser_days_since",
  "N",
  "grass_height",
  "grass_biomass",
  "phytomass_index"
)

# Short labels for plot axes
MGMT_LABELS <- c(
  Grazing_days_since    = "Grazing",
  Fertiliser_days_since = "Fertiliser",
  N                     = "N",
  grass_height          = "Height",
  grass_biomass         = "Biomass",
  phytomass_index       = "PI"
)

# Prediction column suffix per model
PRED_SUFFIX <- c(RF = "rf_predicted", MLP = "mlp_predicted", XGBoost = "xgb_predicted")


# =============================================================================
# SECTION 2 — Libraries
# =============================================================================
# MODEL_COLOURS assigns consistent colours to each model across all six plot
# types.  PRED_COLOURS extends this palette to include each management variable
# and the BASE predictor set, used in the connected-dot plots.
# MGMT_LABELS provides short axis labels for the six management variables.
# MODEL_COLOURS assigns consistent colours to each model across all six plot
# types.  PRED_COLOURS extends this palette to include each management variable
# and the BASE predictor set, used in the connected-dot plots.
# MGMT_LABELS provides short axis labels for the six management variables.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(readr); library(ggplot2); library(tibble)
  library(forcats); library(lubridate)
})

make_dir <- function(...) dir.create(file.path(...), recursive = TRUE, showWarnings = FALSE)

.mae  <- function(y, yh) { ok <- is.finite(y)&is.finite(yh); if(!any(ok)) return(NA_real_); mean(abs(yh[ok]-y[ok])) }
.rmse <- function(y, yh) { ok <- is.finite(y)&is.finite(yh); if(!any(ok)) return(NA_real_); sqrt(mean((yh[ok]-y[ok])^2)) }
.r2   <- function(y, yh) { ok <- is.finite(y)&is.finite(yh); if(sum(ok)<2) return(NA_real_); cc <- suppressWarnings(cor(y[ok],yh[ok])); if(!is.finite(cc)) return(NA_real_); cc^2 }

base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(size = 9),
    legend.position   = "bottom",
    legend.title      = element_text(size = 9),
    legend.text       = element_text(size = 8),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.spacing     = unit(0.8, "lines")
  )

MODEL_COLOURS <- c(RF = "#E41A1C", MLP = "#377EB8", XGBoost = "#4DAF4A")

SEASON_COLOURS <- c(
  Spring = "#4DAF4A", Summer = "#FF7F00",
  Autumn = "#A65628", Winter = "#377EB8", Transitional = "#AAAAAA"
)

metric_ylab <- c(
  mae  = "MAE (\u03bcmol m\u207b\u00b2 s\u207b\u00b9)",
  rmse = "RMSE (\u03bcmol m\u207b\u00b2 s\u207b\u00b9)",
  r2   = "R\u00b2"
)


# =============================================================================
# SECTION 3 — Load BASE metrics from CSVs
# =============================================================================

message("Loading BASE metrics from CSVs ...")
csv_dir <- file.path(GRAPHS_ROOT, "metrics_csv")
if (!dir.exists(csv_dir))
  stop("metrics_csv not found. Run compute_metrics_and_plots.R first.\n", csv_dir)

base_metrics <- map_dfr(TARGET_VARS, function(tv) {
  map_dfr(GAP_SIZES, function(gs) {
    p <- file.path(csv_dir, paste0("gap_metrics_", tv, "_", gs, ".csv"))
    if (!file.exists(p)) return(tibble())
    read_csv(p, show_col_types = FALSE)
  })
}) %>%
  # "unmanaged" in the CSV means BASE predictor set (no management features).
  # Do NOT use management == "managed" — that is a different model run.
  filter(management == "unmanaged", model %in% MODELS) %>%
  mutate(mgmt_var = "BASE",
         gap_size = factor(gap_size, levels = GAP_SIZES))

message("  BASE rows: ", nrow(base_metrics))


# =============================================================================
# SECTION 4 — Compute metrics from management_effect ablation RDS files
# =============================================================================
# compute_mgmt_metrics() mirrors the per-gap metric computation from
# compute_metrics_and_plots.R.  For each gap label in each gap-size category,
# it computes MAE/RMSE/R² for the le30 and gt30 temporal windows.
# Additionally, it classifies each gap label by:
#   recovery_group — fraction of gap rows in the le30 window
#     (High >= 80%, Mid 50–79%, Low < 50%)
#   dominant_season — the season accounting for >50% of gap half-hours
#
# These metadata columns are used by Plots 4 and 5 (recovery interaction and
# seasonal breakdown) to disaggregate the management effect by ecological context.
# The function works identically for all six management variables including PI,
# which is stored in its own subfolder under management_effect/.

message("Computing metrics from management_effect ablation runs ...")

compute_mgmt_metrics <- function(site, model, mgmt_var) {
  
  rds_path <- file.path(RESULTS_ROOT, site, "management_effect",
                        model, mgmt_var, "df_cv_all_predictions.rds")
  if (!file.exists(rds_path)) return(tibble())
  
  preds <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(preds)) return(tibble())
  
  suffix    <- PRED_SUFFIX[[model]]
  flag_all  <- names(preds)[grepl("^(S|M|L|VL)\\d+$", names(preds))]
  has_graz  <- "Grazing_days_since" %in% names(preds)
  has_seas  <- all(c("Spring","Summer","Autumn","Winter") %in% names(preds))
  
  rows <- list()
  for (target in TARGET_VARS) {
    truth_col <- switch(target, NEE="NEE_orig", Reco="Reco_orig", GPP="GPP_orig")
    if (!truth_col %in% names(preds)) next
    
    for (gs in GAP_SIZES) {
      pc <- paste0(target, "_", gs, "_", suffix)
      if (!pc %in% names(preds)) next
      fc_this <- flag_all[str_starts(flag_all, gs)]
      
      for (fc in fc_this) {
        gap_rows <- which(preds[[fc]] %in% c(TRUE, 1L))
        if (!length(gap_rows)) next
        
        y    <- as.numeric(preds[[truth_col]])[gap_rows]
        yhat <- as.numeric(preds[[pc]])[gap_rows]
        sub  <- preds[gap_rows, ]
        
        # Recovery group
        if (has_graz) {
          gd  <- as.numeric(sub$Grazing_days_since)
          pct <- mean(is.finite(gd) & gd <= SPLIT_DAYS, na.rm = TRUE)
        } else pct <- 0
        rec <- case_when(
          pct >= GRAZ_HIGH ~ "High (>=80%)",
          pct >= GRAZ_MID  ~ "Mid (50-79%)",
          TRUE             ~ "Low (<50%)"
        )
        
        # Dominant season
        if (has_seas) {
          props <- c(
            Spring = mean(sub$Spring == 1, na.rm=TRUE),
            Summer = mean(sub$Summer == 1, na.rm=TRUE),
            Autumn = mean(sub$Autumn == 1, na.rm=TRUE),
            Winter = mean(sub$Winter == 1, na.rm=TRUE)
          )
          dom <- if (max(props) > 0.5) names(which.max(props)) else "Transitional"
        } else dom <- "Transitional"
        
        # Split by grazing window (le30 / gt30)
        if (has_graz) {
          gd2  <- as.numeric(preds$Grazing_days_since)[gap_rows]
          s_le <- is.finite(gd2) & gd2 <= SPLIT_DAYS
          s_gt <- is.finite(gd2) & gd2 >  SPLIT_DAYS
        } else {
          s_le <- rep(TRUE,  length(gap_rows))
          s_gt <- rep(FALSE, length(gap_rows))
        }
        
        gap_num <- as.integer(str_extract(fc, "\\d+$"))
        rows[[length(rows)+1]] <- tibble(
          site            = site,
          model           = model,
          mgmt_var        = mgmt_var,
          target          = target,
          gap_size        = gs,
          gap_id          = paste0(gs, gap_num),
          gap_num         = gap_num,
          dominant_season = dom,
          recovery_group  = rec,
          mae_le30  = .mae(y[s_le], yhat[s_le]),
          rmse_le30 = .rmse(y[s_le], yhat[s_le]),
          r2_le30   = .r2(y[s_le], yhat[s_le]),
          mae_gt30  = .mae(y[s_gt], yhat[s_gt]),
          rmse_gt30 = .rmse(y[s_gt], yhat[s_gt]),
          r2_gt30   = .r2(y[s_gt], yhat[s_gt]),
          n_le30    = sum(s_le & is.finite(y) & is.finite(yhat)),
          n_gt30    = sum(s_gt & is.finite(y) & is.finite(yhat))
        )
      }
    }
  }
  if (!length(rows)) return(tibble())
  bind_rows(rows)
}

mgmt_metrics_raw <- map_dfr(SITES, function(site) {
  map_dfr(MODELS, function(model) {
    map_dfr(MGMT_VARS, function(mv) {       # MGMT_VARS now includes "phytomass_index"
      compute_mgmt_metrics(site, model, mv)
    })
  })
}) %>%
  mutate(
    mgmt_var        = factor(mgmt_var, levels = MGMT_VARS),
    gap_size        = factor(gap_size, levels = GAP_SIZES),
    dominant_season = factor(dominant_season,
                             levels = c("Spring","Summer","Autumn","Winter","Transitional")),
    recovery_group  = factor(recovery_group,
                             levels = c("High (>=80%)","Mid (50-79%)","Low (<50%)"))
  )

message("  Management effect metric rows: ", nrow(mgmt_metrics_raw))


# =============================================================================
# SECTION 5 — Reformat BASE metrics for delta computation
# =============================================================================
# The BASE CSV uses the management column value "unmanaged".  Here, the BASE
# metrics are reformatted to match the column structure of mgmt_metrics_raw
# so they can be joined in Section 6.
# IMPORTANT: only rows from the "unmanaged" condition are retained as the BASE
# reference.  The le30/gt30 convention is identical to that in compute_mgmt_metrics
# — both use Grazing_days_since <= SPLIT_DAYS as the le30 definition.
# No column swap is needed.

base_for_join <- base_metrics %>%
  filter(model %in% MODELS) %>%
  rename(gap_size_chr = gap_size) %>%
  mutate(gap_size = factor(gap_size_chr, levels = GAP_SIZES)) %>%
  # Both compute_metrics_one_model and compute_mgmt_metrics define le30/gt30
  # identically (Grazing_days_since <= SPLIT_DAYS).  No swap needed.
  select(site, model, target, gap_size, gap_id,
         mae_le30, rmse_le30, r2_le30,
         mae_gt30, rmse_gt30, r2_gt30,
         n_le30, n_gt30) %>%
  mutate(mgmt_var = "BASE")


# =============================================================================
# SECTION 6 — Aggregate to site × model × management var × gap size and compute deltas
# =============================================================================
# agg_mgmt and agg_base summarise the per-gap-label metrics to one value per
# stratification cell using weighted.mean() with n_le30 or n_gt30 as weights.
# Weighting by observation count ensures that gap labels with more finite pairs
# contribute proportionally more to the aggregate, equivalent to pooling all
# gap rows and computing the metric once.
#
# The delta is then:
#   delta_mae_le30 = agg_mgmt$mae_le30 - agg_base$base_mae_le30
# Negative delta_mae = improvement.  pct_mae_reduction = -100 * delta_mae / base_mae.

agg_mgmt <- mgmt_metrics_raw %>%
  group_by(site, model, mgmt_var, target, gap_size) %>%
  summarise(
    mae_le30  = if (sum(n_le30,  na.rm=TRUE) > 0) weighted.mean(mae_le30,  n_le30,  na.rm=TRUE) else NA_real_,
    rmse_le30 = if (sum(n_le30,  na.rm=TRUE) > 0) weighted.mean(rmse_le30, n_le30,  na.rm=TRUE) else NA_real_,
    r2_le30   = if (sum(n_le30,  na.rm=TRUE) > 0) weighted.mean(r2_le30,   n_le30,  na.rm=TRUE) else NA_real_,
    mae_gt30  = if (sum(n_gt30,  na.rm=TRUE) > 0) weighted.mean(mae_gt30,  n_gt30,  na.rm=TRUE) else NA_real_,
    rmse_gt30 = if (sum(n_gt30,  na.rm=TRUE) > 0) weighted.mean(rmse_gt30, n_gt30,  na.rm=TRUE) else NA_real_,
    r2_gt30   = if (sum(n_gt30,  na.rm=TRUE) > 0) weighted.mean(r2_gt30,   n_gt30,  na.rm=TRUE) else NA_real_,
    .groups   = "drop"
  )

agg_base <- base_for_join %>%
  group_by(site, model, target, gap_size) %>%
  summarise(
    base_mae_le30  = if (sum(n_le30, na.rm=TRUE) > 0) weighted.mean(mae_le30,  n_le30, na.rm=TRUE) else NA_real_,
    base_rmse_le30 = if (sum(n_le30, na.rm=TRUE) > 0) weighted.mean(rmse_le30, n_le30, na.rm=TRUE) else NA_real_,
    base_r2_le30   = if (sum(n_le30, na.rm=TRUE) > 0) weighted.mean(r2_le30,   n_le30, na.rm=TRUE) else NA_real_,
    base_mae_gt30  = if (sum(n_gt30, na.rm=TRUE) > 0) weighted.mean(mae_gt30,  n_gt30, na.rm=TRUE) else NA_real_,
    base_rmse_gt30 = if (sum(n_gt30, na.rm=TRUE) > 0) weighted.mean(rmse_gt30, n_gt30, na.rm=TRUE) else NA_real_,
    base_r2_gt30   = if (sum(n_gt30, na.rm=TRUE) > 0) weighted.mean(r2_gt30,   n_gt30, na.rm=TRUE) else NA_real_,
    .groups = "drop"
  )

delta_metrics <- agg_mgmt %>%
  left_join(agg_base, by = c("site","model","target","gap_size")) %>%
  mutate(
    delta_mae_le30  = mae_le30  - base_mae_le30,
    delta_rmse_le30 = rmse_le30 - base_rmse_le30,
    delta_r2_le30   = r2_le30   - base_r2_le30,
    delta_mae_gt30  = mae_gt30  - base_mae_gt30,
    delta_rmse_gt30 = rmse_gt30 - base_rmse_gt30,
    delta_r2_gt30   = r2_gt30   - base_r2_gt30,
    pct_mae_reduction_le30 = -100 * delta_mae_le30  / base_mae_le30,
    pct_mae_reduction_gt30 = -100 * delta_mae_gt30  / base_mae_gt30,
    mgmt_label = factor(MGMT_LABELS[as.character(mgmt_var)],
                        levels = MGMT_LABELS)
  )

message("Delta metrics computed: ", nrow(delta_metrics), " rows.")


# =============================================================================
# SECTION 7 — Plot 1: Delta-metric barplots — make_delta_barplot()
# =============================================================================
# Grouped bar charts showing ΔMAE, ΔRMSE, or ΔR² for each management variable.
# Bars are grouped by model (fill = model); a horizontal line at y = 0 marks
# no change.  Faceted by site (rows) and gap size (columns) with free y-axes.
# A caption notes the direction of improvement.
# Bars below zero (for MAE/RMSE) or above zero (for R²) indicate improvement.
# Produced separately for the ≤30 d and >30 d windows.

message("\n[1] Delta-metric barplots ...")

make_delta_barplot <- function(tv, metric, window = "le30") {
  delta_col <- paste0("delta_", metric, "_", window)
  win_label <- if (window == "le30") "\u226430 d since grazing" else ">30 d since grazing"
  
  d <- delta_metrics %>%
    filter(target == tv, !is.na(.data[[delta_col]])) %>%
    mutate(
      delta = .data[[delta_col]],
      model = factor(model, levels = MODELS)
    )
  if (!nrow(d)) return(NULL)
  
  ylab <- switch(metric,
                 mae  = "\u0394MAE (\u03bcmol m\u207b\u00b2 s\u207b\u00b9)",
                 rmse = "\u0394RMSE (\u03bcmol m\u207b\u00b2 s\u207b\u00b9)",
                 r2   = "\u0394R\u00b2"
  )
  improve_dir <- if (metric == "r2") "positive = improvement" else "negative = improvement"
  
  ggplot(d, aes(x = mgmt_label, y = delta, fill = model)) +
    geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.9) +
    geom_hline(yintercept = 0, linewidth = 0.6, colour = "black") +
    facet_grid(site ~ gap_size, scales = "free_y") +
    scale_fill_manual(values = MODEL_COLOURS) +
    labs(
      x       = NULL,
      y       = ylab,
      fill    = "Model",
      caption = improve_dir,
      subtitle = paste0(tv, " — ", win_label)
    ) +
    base_theme +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))
}

for (tv in TARGET_VARS) {
  out_dir <- file.path(OUT_ROOT, "delta_barplots", tv)
  make_dir(out_dir)
  for (metric in c("mae","rmse","r2")) {
    for (win in c("le30","gt30")) {
      p <- make_delta_barplot(tv, metric, win)
      if (is.null(p)) next
      ggsave(file.path(out_dir, paste0(toupper(metric), "_", win, ".png")),
             p, width = 14, height = 7, dpi = 220, bg = "white")
    }
  }
  message("  Saved delta barplots: ", tv)
}


# =============================================================================
# SECTION 8 — Plot 2: Connected-dot plots — make_connected_dots()
# =============================================================================
# Shows absolute metric values (not deltas) with the BASE and each management
# variable as separate coloured lines across the three models on the x-axis.
# This layout allows the reader to see both the absolute performance level of
# each predictor set and the shift produced by adding a management variable.
# Faceted by grazing window (rows) and site (columns).

message("\n[2] Connected-dot plots ...")

all_abs <- bind_rows(
  agg_base %>%
    pivot_longer(cols = c(base_mae_le30, base_rmse_le30, base_r2_le30,
                          base_mae_gt30, base_rmse_gt30, base_r2_gt30),
                 names_to  = c("metric","window"),
                 names_pattern = "base_(.+)_(le30|gt30)") %>%
    mutate(predictor_set = "BASE"),
  agg_mgmt %>%
    pivot_longer(cols = c(mae_le30, rmse_le30, r2_le30,
                          mae_gt30, rmse_gt30, r2_gt30),
                 names_to  = c("metric","window"),
                 names_pattern = "(.+)_(le30|gt30)") %>%
    mutate(predictor_set = as.character(mgmt_var))
) %>%
  mutate(
    predictor_set = factor(predictor_set, levels = c("BASE", MGMT_VARS)),
    model         = factor(model, levels = MODELS),
    window_label  = if_else(window == "le30",
                            "\u226430 d since grazing", ">30 d since grazing")
  )

# PI added with a distinct colour (teal)
PRED_COLOURS <- c(
  BASE                  = "#999999",
  Grazing_days_since    = "#E41A1C",
  Fertiliser_days_since = "#FF7F00",
  N                     = "#4DAF4A",
  grass_height          = "#377EB8",
  grass_biomass         = "#984EA3",
  phytomass_index       = "#00BCD4"
)

make_connected_dots <- function(tv, gs, metric) {
  d <- all_abs %>%
    filter(target == tv, gap_size == gs, metric == .env$metric, !is.na(value))
  if (!nrow(d)) return(NULL)
  
  ggplot(d, aes(x = model, y = value,
                colour = predictor_set, group = predictor_set)) +
    geom_line(linewidth = 0.8, alpha = 0.85) +
    geom_point(size = 2.5) +
    facet_grid(window_label ~ site, scales = "free_y") +
    scale_colour_manual(values = PRED_COLOURS,
                        labels = c(MGMT_LABELS, BASE = "BASE"),
                        name   = "Predictor set") +
    labs(x = NULL, y = metric_ylab[[metric]],
         subtitle = paste0(tv, " — Gap size: ", gs)) +
    base_theme +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

for (tv in TARGET_VARS) {
  for (gs in GAP_SIZES) {
    out_dir <- file.path(OUT_ROOT, "connected_dots", tv, gs)
    make_dir(out_dir)
    for (metric in c("mae","rmse","r2")) {
      p <- make_connected_dots(tv, gs, metric)
      if (is.null(p)) next
      ggsave(file.path(out_dir, paste0(toupper(metric), ".png")),
             p, width = 10, height = 7, dpi = 220, bg = "white")
    }
  }
}
message("  Saved connected-dot plots.")


# =============================================================================
# SECTION 9 — Plot 3: Heatmaps of % MAE reduction — make_heatmap()
# =============================================================================
# Two-dimensional summary of the ablation study for one target × gap size:
#   Rows    = management variables (y-axis)
#   Columns = models (x-axis)
#   Cell colour = pct_mae_reduction = -100 * delta_mae / base_mae
# Diverging colour scale: red = degradation, white = no change, green = improvement.
# Percentage values are printed in each cell.  Faceted by site.
# The symmetric colour limit (lim) is rounded up to the nearest 5% so that
# the midpoint is always exactly at 0.

message("\n[3] Heatmaps ...")

make_heatmap <- function(tv, gs, window = "le30") {
  col       <- paste0("pct_mae_reduction_", window)
  win_label <- if (window == "le30") "\u226430 d" else ">30 d"
  
  d <- delta_metrics %>%
    filter(target == tv, gap_size == gs, !is.na(.data[[col]])) %>%
    mutate(
      pct        = .data[[col]],
      model      = factor(model, levels = MODELS),
      mgmt_label = factor(MGMT_LABELS[as.character(mgmt_var)],
                          levels = MGMT_LABELS)
    )
  if (!nrow(d)) return(NULL)
  
  lim <- ceiling(max(abs(d$pct), na.rm = TRUE) / 5) * 5
  
  ggplot(d, aes(x = model, y = mgmt_label, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%+.1f%%", pct)), size = 3) +
    facet_wrap(~ site, nrow = 1) +
    scale_fill_gradient2(
      low  = "#D73027", mid = "white", high = "#1A9850",
      midpoint = 0, limits = c(-lim, lim),
      name = "% MAE reduction\n(positive = improvement)"
    ) +
    labs(x = NULL, y = NULL,
         subtitle = paste0(tv, " — ", gs, " gaps — ", win_label)) +
    base_theme +
    theme(legend.position = "right")
}

for (tv in TARGET_VARS) {
  out_dir <- file.path(OUT_ROOT, "heatmaps", tv)
  make_dir(out_dir)
  for (gs in GAP_SIZES) {
    for (win in c("le30","gt30")) {
      p <- make_heatmap(tv, gs, win)
      if (is.null(p)) next
      ggsave(file.path(out_dir, paste0("pct_MAE_", gs, "_", win, ".png")),
             p, width = 9, height = 5, dpi = 220, bg = "white")
    }
  }
}
message("  Saved heatmaps.")


# =============================================================================
# SECTION 10 — Plot 4: Recovery-phase interaction — make_recovery_plot()
# =============================================================================
# Tests whether the management effect is stronger during active regrowth
# (early post-grazing) than during later recovery phases.  For each gap label,
# the recovery_group (High/Mid/Low) was assigned in Section 4.  The delta metric
# is averaged within each recovery group and plotted as connected lines across
# groups, with one line per management variable.  Faceted by model and site.
# A dashed reference line at y = 0 marks no effect.

message("\n[4] Recovery-phase interaction plots ...")

recovery_delta <- mgmt_metrics_raw %>%
  left_join(
    base_for_join %>%
      select(site, model, target, gap_size, gap_id,
             base_mae_le30 = mae_le30, base_mae_gt30 = mae_gt30,
             base_r2_le30  = r2_le30,  base_r2_gt30  = r2_gt30),
    by = c("site","model","target","gap_size","gap_id")
  ) %>%
  mutate(
    delta_mae_le30 = mae_le30 - base_mae_le30,
    delta_mae_gt30 = mae_gt30 - base_mae_gt30,
    delta_r2_le30  = r2_le30  - base_r2_le30,
    delta_r2_gt30  = r2_gt30  - base_r2_gt30
  )

make_recovery_plot <- function(tv, metric) {
  delta_col <- paste0("delta_", metric, "_le30")
  
  d <- recovery_delta %>%
    filter(target == tv, !is.na(.data[[delta_col]]), !is.na(recovery_group)) %>%
    group_by(site, model, mgmt_var, recovery_group) %>%
    summarise(mean_delta = mean(.data[[delta_col]], na.rm = TRUE), .groups = "drop") %>%
    mutate(
      mgmt_label     = factor(MGMT_LABELS[as.character(mgmt_var)], levels = MGMT_LABELS),
      recovery_group = factor(recovery_group,
                              levels = c("High (>=80%)","Mid (50-79%)","Low (<50%)")),
      model = factor(model, levels = MODELS)
    )
  if (!nrow(d)) return(NULL)
  
  ylab <- switch(metric,
                 mae = "\u0394MAE (negative = improvement)",
                 r2  = "\u0394R\u00b2 (positive = improvement)"
  )
  
  ggplot(d, aes(x = recovery_group, y = mean_delta,
                colour = mgmt_label, group = mgmt_label)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    facet_grid(model ~ site) +
    labs(x = "Grazing recovery phase", y = ylab,
         colour = "Management variable",
         subtitle = paste0(tv, " — \u226430 d window")) +
    base_theme +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
}

for (tv in TARGET_VARS) {
  out_dir <- file.path(OUT_ROOT, "recovery_interaction", tv)
  make_dir(out_dir)
  for (metric in c("mae","r2")) {
    p <- make_recovery_plot(tv, metric)
    if (is.null(p)) next
    ggsave(file.path(out_dir, paste0(toupper(metric), ".png")),
           p, width = 10, height = 8, dpi = 220, bg = "white")
  }
}
message("  Saved recovery-interaction plots.")


# =============================================================================
# SECTION 11 — Plot 5: Seasonal breakdown — make_seasonal_plot()
# =============================================================================
# Tests whether management variables are more or less useful depending on the
# season in which the gap falls.  Boxplots of ΔMAE for each management variable,
# with boxes coloured by the dominant season of each gap label.  Faceted by
# model and site.  A horizontal line at y = 0 marks no effect.
# Uses the ≤30 d grazing window only (where management effects are strongest).

message("\n[5] Seasonal breakdown ...")

seasonal_delta <- mgmt_metrics_raw %>%
  left_join(
    base_for_join %>%
      select(site, model, target, gap_size, gap_id,
             base_mae_le30 = mae_le30),
    by = c("site","model","target","gap_size","gap_id")
  ) %>%
  mutate(delta_mae = mae_le30 - base_mae_le30)

make_seasonal_plot <- function(tv) {
  d <- seasonal_delta %>%
    filter(target == tv, !is.na(delta_mae), !is.na(dominant_season)) %>%
    mutate(
      mgmt_label = factor(MGMT_LABELS[as.character(mgmt_var)], levels = MGMT_LABELS),
      model      = factor(model, levels = MODELS)
    )
  if (!nrow(d)) return(NULL)
  
  ggplot(d, aes(x = mgmt_label, y = delta_mae, fill = dominant_season)) +
    geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5,
                 position = position_dodge(0.8), width = 0.65) +
    geom_hline(yintercept = 0, linewidth = 0.5, colour = "black") +
    facet_grid(model ~ site) +
    scale_fill_manual(values = SEASON_COLOURS, name = "Dominant season") +
    labs(x = NULL,
         y = "\u0394MAE (\u03bcmol m\u207b\u00b2 s\u207b\u00b9) | negative = improvement",
         subtitle = paste0(tv, " — \u226430 d window")) +
    base_theme +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))
}

for (tv in TARGET_VARS) {
  out_dir <- file.path(OUT_ROOT, "seasonal_breakdown", tv)
  make_dir(out_dir)
  p <- make_seasonal_plot(tv)
  if (!is.null(p))
    ggsave(file.path(out_dir, "delta_MAE_by_season.png"),
           p, width = 12, height = 9, dpi = 220, bg = "white")
}
message("  Saved seasonal breakdown plots.")


# =============================================================================
# SECTION 12 — Plot 6: Ranking stability bump charts — make_bump_chart()
# =============================================================================
# A bump chart shows the rank order of management variables (y-axis, reversed
# so rank 1 = best = top) across the three models (x-axis).  Each management
# variable is a coloured line; rank numbers are printed at each model position.
# Parallel lines = consistent ranking across models.  Crossing lines = the
# ranking depends on the model choice, suggesting model-specific interactions
# with certain management variables.
# Produced for MAE and R² separately, and for both grazing windows.

message("\n[6] Ranking stability bump charts ...")

make_bump_chart <- function(tv, gs, metric, window = "le30") {
  delta_col <- paste0("delta_", metric, "_", window)
  
  ranks <- delta_metrics %>%
    filter(target == tv, gap_size == gs, !is.na(.data[[delta_col]])) %>%
    group_by(site, model) %>%
    mutate(
      rank = if (metric == "r2")
        rank(-(.data[[delta_col]]), ties.method = "min")
      else
        rank(.data[[delta_col]], ties.method = "min"),
      mgmt_label = factor(MGMT_LABELS[as.character(mgmt_var)],
                          levels = MGMT_LABELS),
      model = factor(model, levels = MODELS)
    ) %>%
    ungroup()
  
  if (!nrow(ranks)) return(NULL)
  
  n_vars    <- length(MGMT_VARS)
  win_label <- if (window == "le30") "\u226430 d" else ">30 d"
  
  ggplot(ranks, aes(x = model, y = rank,
                    colour = mgmt_label, group = mgmt_label)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 3) +
    geom_text(aes(label = rank), size = 2.5, colour = "white", fontface = "bold") +
    scale_y_reverse(breaks = 1:n_vars,
                    labels = paste0("Rank ", 1:n_vars)) +
    facet_wrap(~ site, nrow = 1) +
    labs(x = NULL, y = "Rank (1 = best improvement)",
         colour = "Management variable",
         subtitle = paste0(tv, " — ", gs, " gaps — ",
                           toupper(metric), " — ", win_label)) +
    base_theme +
    theme(panel.grid.major.x = element_blank())
}

for (tv in TARGET_VARS) {
  out_dir <- file.path(OUT_ROOT, "ranking_stability", tv)
  make_dir(out_dir)
  for (gs in GAP_SIZES) {
    for (metric in c("mae","r2")) {
      p <- make_bump_chart(tv, gs, metric)
      if (is.null(p)) next
      ggsave(file.path(out_dir, paste0(metric, "_", gs, ".png")),
             p, width = 8, height = 5, dpi = 220, bg = "white")
    }
  }
}
message("  Saved ranking stability bump charts.")


# =============================================================================
# SECTION 13 — Summary
# =============================================================================

message("\n", strrep("=", 60))
message("ALL MANAGEMENT EFFECT PLOTS COMPLETE")
message(strrep("=", 60))
message("Output root: ", OUT_ROOT)
message("  1. delta_barplots/        — ΔMAE/ΔRMSE/ΔR² per variable, fill=model")
message("  2. connected_dots/        — absolute metric lines by predictor set")
message("  3. heatmaps/              — % MAE reduction, rows=var, cols=model")
message("  4. recovery_interaction/  — delta by early vs late recovery phase")
message("  5. seasonal_breakdown/    — delta boxplots by dominant season")
message("  6. ranking_stability/     — bump charts of variable rank across models")