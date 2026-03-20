# =============================================================================
# compute_metrics_and_plots.R — Gap-Filling Metrics and Visualisation
# =============================================================================
#
# PURPOSE
#   Reads df_cv_all_predictions.rds from every model result folder, computes
#   per-gap-label quality metrics (MAE, RMSE, R²), and produces two families
#   of standard evaluation plots.  The metrics are also exported as CSV files
#   for use by the management effect and variable importance plotting scripts.
#
# METRICS COMPUTED
#   For every site × model × management × target × gap size × gap label:
#     mae_all, rmse_all, r2_all   — computed over all gap rows
#     mae_le30, rmse_le30, r2_le30 — restricted to gap rows where
#                                    Grazing_days_since <= SPLIT_DAYS
#     mae_gt30, rmse_gt30, r2_gt30 — restricted to gap rows where
#                                    Grazing_days_since >  SPLIT_DAYS
#   The grazing-window split is ecologically motivated: early post-grazing
#   periods (<=30 d) have rapidly changing canopy structure and carbon balance,
#   making gap-filling systematically harder.  Reporting metrics separately
#   for both windows reveals whether management-informed models have their
#   greatest advantage in the early recovery phase.
#
# PLOTS PRODUCED
#   (1) OVERALL METRICS  →  graphs/overall_metrics/{TARGET}/
#       Line-profile plots: x = gap size (S → VL), y = metric value,
#       coloured by model.  Rows = ≤30 d vs >30 d windows.
#       Columns = JC1 Unmanaged | JC1 Managed | JC2 Unmanaged | JC2 Managed.
#       Saved as PNG (no title) + PDF (with title), one per metric.
#
#   (2) SITE COMPARISON  →  graphs/site_comparison/{TARGET}/{GAP_SIZE}/
#       Connected-dot plots: x = model (ordered best → worst), y = metric.
#       Each grey line is one artificial gap label; coloured dots are
#       individual gap values.  The black diamond line is the weighted-mean
#       aggregate across all gap labels (weighted by n_le30 or n_gt30).
#       Saved as PNG + PDF, one per metric × gap size × target.
#
# INPUTS
#   results/{SITE}/{MANAGEMENT}/{MODEL}/df_cv_all_predictions.rds
#
# OUTPUTS
#   graphs/metrics_csv/gap_metrics_{TARGET}_{GAP_SIZE}.csv  — raw per-gap metrics
#   graphs/overall_metrics/{TARGET}/                        — line-profile plots
#   graphs/site_comparison/{TARGET}/{GAP_SIZE}/             — connected-dot plots
#
# SECTIONS
#   1  User settings       — site, split threshold, directories, target vars
#   2  Package checks      — verify and load required libraries; define helpers
#   3  Load and compute    — read RDS files, compute per-gap metrics, save CSVs
#   4  Plot theme          — shared ggplot2 theme, colour palette, axis labels
#   5  Overall metrics     — make_overall_plot() and save loop
#   6  Site comparison     — make_comparison_plot() and save loop
#   8  Summary message
#
# HOW TO RUN
#   source("metrics/compute_metrics_and_plots.R")
#   Run this script BEFORE management_effect_graphs.R and
#   variable_importance_graphs.R, which read the CSVs it produces.
#
# =============================================================================


# =============================================================================
# SECTION 1 — User settings
# =============================================================================
# Edit the variables in this block to adapt the script to a different site
# or analysis configuration.  All other sections are automatic.

# --- Site to analyse in overall metrics plots --------------------------------
SITE <- "JC1"    # change to "JC2" for the second site

# --- Grazing-window split threshold (days) -----------------------------------
# Gap half-hours with Grazing_days_since <= SPLIT_DAYS are classified as early
# recovery (le30); those above the threshold are late recovery (gt30).
SPLIT_DAYS <- 30    # 30 days is the standard post-grazing split

# --- Root directories --------------------------------------------------------
RESULTS_ROOT <- here::here("results")   # where df_cv_all_predictions.rds files live
GRAPHS_ROOT  <- here::here("graphs")    # where all output plots will be saved

# --- Target variables to process ---------------------------------------------
TARGET_VARS <- c("NEE", "Reco", "GPP")


# =============================================================================
# SECTION 2 — Package checks and helpers
# =============================================================================
# All required packages are checked for availability before any computation
# begins.  The model catalogue and path-building functions defined here are
# shared with graph_colour_individual_gaps.R; both scripts must use the same
# MODEL_CATALOG so that model labels and management assignments are consistent.

pkgs <- c("here","dplyr","tibble","tidyr","purrr","stringr","readr",
          "ggplot2","patchwork","grDevices")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tibble); library(tidyr)
  library(purrr); library(stringr); library(readr); library(ggplot2)
  library(patchwork)
})

# ---- Metric helpers ---------------------------------------------------------
# All three functions guard against NA/Inf in both observed and predicted
# values.  Only finite pairs are used; NA is returned when no finite pairs
# exist rather than propagating NaN through downstream summaries.
.mae  <- function(y, yh) { ok <- is.finite(y)&is.finite(yh); if(!any(ok)) return(NA_real_); mean(abs(yh[ok]-y[ok])) }
.rmse <- function(y, yh) { ok <- is.finite(y)&is.finite(yh); if(!any(ok)) return(NA_real_); sqrt(mean((yh[ok]-y[ok])^2)) }
.r2   <- function(y, yh) { ok <- is.finite(y)&is.finite(yh); if(sum(ok)<2) return(NA_real_); cc <- suppressWarnings(cor(y[ok],yh[ok])); if(!is.finite(cc)) return(NA_real_); cc^2 }

# ---- Model catalogue --------------------------------------------------------
# Maps each model to its result folder name (model_key), display label
# (model_label), and management condition (management).  miniRECgap and MDS
# live outside the managed/unmanaged folder hierarchy (management = NA) and
# their RDS files are located directly under results/{SITE}/{MODEL_NAME}/.
# This table drives both the metric computation loop and the path builder.
MODEL_CATALOG <- tribble(
  ~model_key,  ~model_label,  ~management,
  "rf",        "RF",          "managed",
  "mlp",       "MLP",         "managed",
  "xgb",       "XGBoost",     "managed",
  "rf",        "RF",          "unmanaged",
  "mlp",       "MLP",         "unmanaged",
  "xgb",       "XGBoost",     "unmanaged",
  "minirec",   "miniRECgap",  NA_character_,
  "mds",       "MDS",         NA_character_
)

# ---- Prediction-column lookup — pred_col_name() ----------------------------
# Constructs the prediction column name used by each model in
# df_cv_all_predictions.rds, e.g. NEE_M_rf_predicted.
pred_col_name <- function(target, gap_size, model_key) {
  suffix <- switch(model_key,
                   rf      = "rf_predicted",
                   mlp     = "mlp_predicted",
                   xgb     = "xgb_predicted",
                   minirec = "minirec_predicted",
                   mds     = "mds_predicted",
                   stop("Unknown model key: ", model_key))
  paste0(target, "_", gap_size, "_", suffix)
}

# ---- Path builder — rds_path_for_model() ------------------------------------
# Constructs the full path to df_cv_all_predictions.rds for a given
# site × model × management combination.  miniRECgap and MDS use a flat
# folder structure (results/{SITE}/{MODEL}/); the three ML models use a
# nested structure (results/{SITE}/{management}/{MODEL}/).
rds_path_for_model <- function(site, model_key, management) {
  if (model_key %in% c("minirec","mds")) {
    folder <- switch(model_key, minirec="miniRECgap", mds="MDS")
    file.path(RESULTS_ROOT, site, folder, "df_cv_all_predictions.rds")
  } else {
    folder <- switch(model_key, rf="RF", mlp="MLP", xgb="XGBoost")
    file.path(RESULTS_ROOT, site, management, folder, "df_cv_all_predictions.rds")
  }
}

# ---- Safe directory creator -------------------------------------------------
make_dir <- function(...) dir.create(file.path(...), recursive=TRUE, showWarnings=FALSE)


# =============================================================================
# SECTION 3 — Load all predictions and compute metrics
# =============================================================================
# Iterates over the full MODEL_CATALOG × SITES_ALL grid, reading each model's
# df_cv_all_predictions.rds and computing per-gap-label quality metrics.
# Three metric windows are computed for each gap label:
#   (a) All gap rows (mae_all / rmse_all / r2_all)
#   (b) Gap rows where Grazing_days_since ≤ SPLIT_DAYS  (le30 metrics)
#   (c) Gap rows where Grazing_days_since >  SPLIT_DAYS  (gt30 metrics)
# The window split requires Grazing_days_since in the RDS; if absent, all
# rows are treated as "le30" and "gt30" returns all-NA metrics.
# n_le30 and n_gt30 record the number of finite observation pairs in each
# window, used as weights in the aggregate summaries.
# Results are stored in the flat tibble `all_gap_metrics` and immediately
# exported to CSV for use by the management effect and VI plotting scripts.

message("Loading predictions and computing gap-level metrics ...")

compute_metrics_one_model <- function(site, model_key, management) {
  
  path <- rds_path_for_model(site, model_key, management)
  if (!file.exists(path)) {
    message("  Not found (skipping): ", path); return(tibble())
  }
  
  preds <- tryCatch(readRDS(path), error = function(e) {
    message("  Failed to load ", path, ": ", conditionMessage(e)); return(NULL)
  })
  if (is.null(preds)) return(tibble())
  
  if (!"NEE_orig" %in% names(preds)) {
    message("  NEE_orig missing in ", path); return(tibble())
  }
  
  has_grazing <- "Grazing_days_since" %in% names(preds)
  if (!has_grazing)
    message("  Grazing_days_since absent in ", path,
            " — all rows treated as one group.")
  
  model_label <- MODEL_CATALOG %>%
    filter(
      model_key == .env$model_key,
      (is.na(.env$management) & is.na(management)) |
        (!is.na(.env$management) & !is.na(management) & management == .env$management)
    ) %>%
    pull(model_label) %>% first()
  if (is.na(model_label)) model_label <- model_key
  
  rows <- list()
  
  for (target in TARGET_VARS) {
    for (gap_size in c("S","M","L","VL")) {
      
      # Locate gap membership flags for this gap-size category
      flag_cols <- names(preds)[grepl(paste0("^", gap_size, "\\d+$"), names(preds))]
      if (!length(flag_cols)) next
      
      pc <- pred_col_name(target, gap_size, model_key)
      if (!pc %in% names(preds)) next
      
      # Ground-truth column selection:
      # NEE uses the directly measured NEE_orig.
      # Reco_orig and GPP_orig are derived by flux partitioning (Lloyd-Taylor /
      # Thornley) applied to NEE_orig per regrowth period by the model scripts.
      # They are not independently measured and their accuracy depends on both
      # the quality of NEE_orig and the fit of the partitioning parameters.
      truth_col <- switch(target,
                          NEE  = "NEE_orig",
                          Reco = "Reco_orig",
                          GPP  = "GPP_orig",
                          "NEE_orig"   # safe fallback
      )
      if (!truth_col %in% names(preds)) {
        message("  Ground-truth column '", truth_col,
                "' missing for target=", target, " in ", path,
                " — run patch_reco_gpp.R first. Skipping.")
        next
      }
      
      for (fc in flag_cols) {
        gap_rows <- which(preds[[fc]] %in% c(TRUE, 1))
        if (!length(gap_rows)) next
        
        y    <- as.numeric(preds[[truth_col]])[gap_rows]
        yhat <- as.numeric(preds[[pc]])[gap_rows]
        
        if (has_grazing) {
          gd <- as.numeric(preds$Grazing_days_since)[gap_rows]
          sel_le <- is.finite(gd) & gd <= SPLIT_DAYS
          sel_gt <- is.finite(gd) & gd >  SPLIT_DAYS
        } else {
          sel_le <- rep(TRUE,  length(gap_rows))
          sel_gt <- rep(FALSE, length(gap_rows))
        }
        
        gap_num <- as.integer(str_extract(fc, "\\d+$"))
        rows[[length(rows)+1]] <- tibble(
          site       = site,
          management = coalesce(management, "unmanaged"),
          model      = model_label,
          model_key  = model_key,
          target     = target,
          gap_size   = gap_size,
          gap_id     = paste0(gap_size, gap_num),
          gap_num    = gap_num,
          mae_all    = .mae(y, yhat),
          rmse_all   = .rmse(y, yhat),
          r2_all     = .r2(y, yhat),
          n_le30     = sum(sel_le & is.finite(y) & is.finite(yhat)),
          n_gt30     = sum(sel_gt & is.finite(y) & is.finite(yhat)),
          mae_le30   = .mae(y[sel_le], yhat[sel_le]),
          rmse_le30  = .rmse(y[sel_le], yhat[sel_le]),
          r2_le30    = .r2(y[sel_le], yhat[sel_le]),
          mae_gt30   = .mae(y[sel_gt], yhat[sel_gt]),
          rmse_gt30  = .rmse(y[sel_gt], yhat[sel_gt]),
          r2_gt30    = .r2(y[sel_gt], yhat[sel_gt])
        )
      }
    }
  }
  if (!length(rows)) return(tibble())
  bind_rows(rows)
}

# Process all combinations
SITES_ALL <- c("JC1", "JC2")

all_gap_metrics <- pmap_dfr(
  MODEL_CATALOG,
  function(model_key, model_label, management) {
    map_dfr(SITES_ALL, function(site) {
      compute_metrics_one_model(site, model_key, management)
    })
  }
) %>%
  mutate(
    model     = factor(model, levels = c("RF","MLP","XGBoost","miniRECgap","MDS")),
    gap_size  = factor(gap_size, levels = c("S","M","L","VL")),
    management = factor(management, levels = c("unmanaged","managed"))
  )

message("Gap-level metrics computed: ", nrow(all_gap_metrics), " rows.")

# ---- Save CSVs — one file per target × gap size ----------------------------
# These CSVs are the primary output of this script and are read by
# management_effect_graphs.R and graph_colour_individual_gaps.R.
# Each file contains the full per-gap-label metric table including the
# le30/gt30 window split, n_le30, n_gt30, site, model, and management columns.
csv_dir <- file.path(GRAPHS_ROOT, "metrics_csv")
make_dir(csv_dir)
for (tv in TARGET_VARS) {
  for (gs in c("S","M","L","VL")) {
    d <- filter(all_gap_metrics, target == tv, gap_size == gs)
    write_csv(d, file.path(csv_dir, paste0("gap_metrics_", tv, "_", gs, ".csv")))
  }
}
message("CSVs saved to: ", csv_dir)

# ---- Drop MDS from long gaps ------------------------------------------------
# MDS requires a sufficient density of similar meteorological conditions within
# its search window.  For very long gaps (VL = ~30 days, L = ~14 days), the
# search window can exhaust available analogues, producing unreliable fills.
# MDS results for L and VL gap sizes are therefore excluded from all plots.
all_gap_metrics <- all_gap_metrics %>%
  filter(!(model == "MDS" & gap_size %in% c("L", "VL")))

# ---- overall_metrics: simple mean across gap instances ----------------------
# Used by make_overall_plot() (Section 5) to draw the model line profiles.
# Simple unweighted mean is used here because each gap label is treated as one
# observation of model performance at a given gap size.
overall_metrics <- all_gap_metrics %>%
  group_by(site, management, model, model_key, target, gap_size) %>%
  summarise(
    mae_le30  = mean(mae_le30,  na.rm=TRUE),
    rmse_le30 = mean(rmse_le30, na.rm=TRUE),
    r2_le30   = mean(r2_le30,   na.rm=TRUE),
    mae_gt30  = mean(mae_gt30,  na.rm=TRUE),
    rmse_gt30 = mean(rmse_gt30, na.rm=TRUE),
    r2_gt30   = mean(r2_gt30,   na.rm=TRUE),
    .groups = "drop"
  )

# ---- agg_metrics: weighted mean across gap instances ------------------------
# Used by make_comparison_plot() (Section 6) as the black aggregate reference
# line.  Weighted by n_le30 or n_gt30 so that longer or data-richer gaps
# contribute proportionally more to the aggregate than short data-sparse ones.
# This produces the same metric one would obtain by pooling all gap rows into
# a single set and computing MAE/RMSE/R² once.
agg_metrics <- all_gap_metrics %>%
  group_by(site, management, model, model_key, target, gap_size) %>%
  summarise(
    agg_mae_le30  = weighted.mean(mae_le30,  n_le30,  na.rm=TRUE),
    agg_rmse_le30 = weighted.mean(rmse_le30, n_le30,  na.rm=TRUE),
    agg_r2_le30   = weighted.mean(r2_le30,   n_le30,  na.rm=TRUE),
    agg_mae_gt30  = weighted.mean(mae_gt30,  n_gt30,  na.rm=TRUE),
    agg_rmse_gt30 = weighted.mean(rmse_gt30, n_gt30,  na.rm=TRUE),
    agg_r2_gt30   = weighted.mean(r2_gt30,   n_gt30,  na.rm=TRUE),
    .groups = "drop"
  )


# =============================================================================
# SECTION 4 — Plot theme and shared helpers
# =============================================================================
# A single ggplot2 theme (base_theme) is applied to all figures in this script
# for visual consistency.  MODEL_COLOURS assigns the same colour to each model
# across all plot types; MDS is grey because it is included only as a reference.
# window_labels provides human-readable row facet labels for the ≤30 d / >30 d
# split.  metric_ylab provides Unicode axis labels with proper units.

base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor    = element_blank(),
    strip.text          = element_text(size = 11),
    legend.position     = "bottom",
    legend.title        = element_text(size = 9),
    legend.text         = element_text(size = 8),
    plot.background     = element_rect(fill = "white", colour = NA),
    panel.spacing       = unit(1, "lines")
  )

# Colour palette (colour-blind friendly)
MODEL_COLOURS <- c(
  RF          = "#E41A1C",
  MLP         = "#377EB8",
  XGBoost     = "#4DAF4A",
  miniRECgap  = "#FF7F00",
  MDS         = "#999999"   # grey = empty placeholder
)

window_labels <- c(
  le30 = paste0("\u226430 d since grazing"),
  gt30 = paste0(">30 d since grazing")
)

metric_ylab <- c(
  mae  = "MAE (\u03bcmol m\u207b\u00b2 s\u207b\u00b9)",
  rmse = "RMSE (\u03bcmol m\u207b\u00b2 s\u207b\u00b9)",
  r2   = "R\u00b2"
)



# =============================================================================
# SECTION 5 — Overall metrics plots
# =============================================================================
# make_overall_plot() draws a 4-column × 2-row faceted line chart:
#   Columns: JC1 Unmanaged | JC1 Managed | JC2 Unmanaged | JC2 Managed
#   Rows:    ≤30 d since grazing  |  >30 d since grazing
#   x-axis:  gap size (S → M → L → VL)
#   y-axis:  the selected metric (MAE, RMSE, or R²)
#   Colour:  model identity
# scales = "free_y" allows each row to share one y-axis range across all four
# columns, enabling direct visual comparison between sites and management
# conditions without one panel distorting the others.
# MDS appears only in the Unmanaged columns (no managed run exists).
# Plots are saved as PNG (no title, for embedding in figures) and as a
# multi-page PDF with titles (for internal review).

message("\nBuilding OVERALL METRICS plots (4-column: JC1/JC2 × managed/unmanaged) ...")

COL_LABEL_LEVELS <- c("JC1 Unmanaged","JC1 Managed","JC2 Unmanaged","JC2 Managed")

make_overall_plot <- function(tv, metric, add_title = FALSE, title_str = "") {
  
  d <- overall_metrics %>%
    filter(target == tv) %>%
    mutate(
      col_label = factor(paste(site, str_to_title(management)),
                         levels = COL_LABEL_LEVELS)
    ) %>%
    select(col_label, model, gap_size,
           value_le30 = paste0(metric, "_le30"),
           value_gt30 = paste0(metric, "_gt30")) %>%
    pivot_longer(c(value_le30, value_gt30),
                 names_to = "window", values_to = "value",
                 names_pattern = "value_(le30|gt30)") %>%
    mutate(window = factor(window, levels = c("le30","gt30"),
                           labels = window_labels)) %>%
    filter(!is.na(value))
  
  if (!nrow(d)) return(NULL)
  
  p <- ggplot(d, aes(x = gap_size, y = value,
                     colour = model, group = model)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    facet_grid(window ~ col_label, scales = "free_y") +
    scale_colour_manual(values = MODEL_COLOURS, drop = FALSE) +
    labs(x = "Gap size", y = metric_ylab[[metric]],
         colour = "Model",
         title  = if (add_title) title_str else NULL) +
    base_theme +
    theme(
      legend.position = "bottom",
      strip.text.y    = element_text(angle = -90, size = 9)
    )
  
  y_vals <- d$value[is.finite(d$value)]           # ← d is the data frame in this function
  if (length(y_vals)) {
    if (metric == "r2") {
      y_lo <- floor(min(y_vals) * 10) / 10    # round down to nearest 0.1
      y_hi <- ceiling(max(y_vals) * 10) / 10  # round up to nearest 0.1
      p <- p + coord_cartesian(ylim = c(y_lo, y_hi)) +
        scale_y_continuous(breaks = seq(y_lo, y_hi, by = 0.1))
    } else {
      p <- p + scale_y_continuous(
        breaks = seq(floor(min(y_vals)), ceiling(max(y_vals)), by = 1)
      )
    }
  }
  
  p
}

for (tv in TARGET_VARS) {
  out_dir <- file.path(GRAPHS_ROOT, "overall_metrics", tv)
  make_dir(out_dir)
  
  pdf_plots <- list()
  for (metric in c("mae","rmse","r2")) {
    
    p_png <- make_overall_plot(tv, metric, add_title = FALSE)
    if (is.null(p_png)) next
    out_png <- file.path(out_dir, paste0(toupper(metric), "_overall.png"))
    ggsave(out_png, p_png, width = 14, height = 7, dpi = 220, bg = "white")
    message("  Saved: ", out_png)
    
    pdf_plots[[metric]] <- make_overall_plot(
      tv, metric, add_title = TRUE,
      title_str = paste0(tv, " — ", toupper(metric),
                         " by gap size | JC1 & JC2, managed & unmanaged"))
  }
  
  out_pdf <- file.path(out_dir, paste0("overall_metrics_", tv, ".pdf"))
  pdf(out_pdf, width = 14, height = 7, onefile = TRUE)
  for (p in pdf_plots) if (!is.null(p)) print(p)
  dev.off()
  message("  Saved PDF: ", out_pdf)
}


# =============================================================================
# SECTION 6 — Site comparison (connected-dot plots)
# =============================================================================
# make_comparison_plot() draws the individual-gap connected-dot plot with the
# following design:
#   - Each point = one artificial gap label's metric value.
#   - Grey lines connect the same gap label across models, showing which
#     gaps are systematically easy or hard for all models.
#   - Coloured dots show per-model per-gap values.
#   - The black diamond line is the weighted-aggregate reference.
#   - Models on the x-axis are ordered best → worst using JC1 Unmanaged
#     as the reference column, so the best-performing model is leftmost.
# model_order_global() computes this ordering from agg_metrics.
# mods_for() returns the correct model set for each management condition
# and gap size (MDS excluded from L and VL, managed columns exclude miniRECgap
# and MDS which have no managed run).
# Both PNG (no title) and PDF (with title) are saved for each combination.

message("\nBuilding SITE COMPARISON / individual gap plots ...")

# Helper: compute global model order from JC1 unmanaged aggregate
model_order_global <- function(tv, gap_cat, metric) {
  higher_better <- metric == "r2"
  ref_col <- paste0("agg_", metric, "_le30")
  agg_metrics %>%
    filter(target == tv, gap_size == gap_cat,
           site == "JC1", management == "unmanaged",
           !is.na(.data[[ref_col]])) %>%
    arrange(if (higher_better) desc(.data[[ref_col]]) else .data[[ref_col]]) %>%
    pull(model) %>% as.character()
}

# Models shown per management × gap_size
mods_for <- function(mgmt, gap_cat) {
  if (mgmt == "managed")          return(c("RF","MLP","XGBoost"))
  if (gap_cat %in% c("VL","L"))  return(c("RF","MLP","XGBoost","miniRECgap"))
  return(c("RF","MLP","XGBoost","miniRECgap","MDS"))
}

make_comparison_plot <- function(tv, gap_cat, metric,
                                 add_title = FALSE, title_str = "") {
  value_le <- paste0(metric, "_le30")
  value_gt <- paste0(metric, "_gt30")
  agg_le   <- paste0("agg_", metric, "_le30")
  agg_gt   <- paste0("agg_", metric, "_gt30")
  
  # Determine global model order from JC1 unmanaged reference
  ord <- model_order_global(tv, gap_cat, metric)
  # Append any remaining models not in reference (e.g. from managed only)
  all_mods <- unique(c(
    mods_for("unmanaged", gap_cat),
    mods_for("managed",   gap_cat)
  ))
  ord_full <- c(ord, setdiff(all_mods, ord))
  
  # Build individual gap data (long: one row per gap × model × window)
  col_specs <- list(
    list(site="JC1", mgmt="unmanaged"),
    list(site="JC1", mgmt="managed"),
    list(site="JC2", mgmt="unmanaged"),
    list(site="JC2", mgmt="managed")
  )
  
  d_indiv <- map_dfr(col_specs, function(spec) {
    site <- spec$site; mgmt <- spec$mgmt
    mods <- mods_for(mgmt, gap_cat)
    all_gap_metrics %>%
      filter(target==tv, gap_size==gap_cat,
             site==.env$site, management==mgmt,
             model %in% mods) %>%
      select(model, gap_id,
             le30 = all_of(value_le),
             gt30 = all_of(value_gt)) %>%
      pivot_longer(c(le30, gt30), names_to="window", values_to="value") %>%
      mutate(col_label = paste(site, str_to_title(mgmt)))
  }) %>%
    filter(!is.na(value)) %>%
    mutate(
      window    = factor(window, levels=c("le30","gt30"), labels=window_labels),
      col_label = factor(col_label, levels=COL_LABEL_LEVELS),
      model     = factor(model, levels=ord_full)
    )
  
  # Build aggregate reference data
  d_agg <- map_dfr(col_specs, function(spec) {
    site <- spec$site; mgmt <- spec$mgmt
    mods <- mods_for(mgmt, gap_cat)
    agg_metrics %>%
      filter(target==tv, gap_size==gap_cat,
             site==.env$site, management==mgmt,
             model %in% mods) %>%
      select(model,
             le30 = all_of(agg_le),
             gt30 = all_of(agg_gt)) %>%
      pivot_longer(c(le30, gt30), names_to="window", values_to="agg_value") %>%
      mutate(col_label = paste(site, str_to_title(mgmt)))
  }) %>%
    filter(!is.na(agg_value)) %>%
    mutate(
      window    = factor(window, levels=c("le30","gt30"), labels=window_labels),
      col_label = factor(col_label, levels=COL_LABEL_LEVELS),
      model     = factor(model, levels=ord_full)
    )
  
  if (!nrow(d_indiv)) return(NULL)
  
  p <- ggplot(d_indiv, aes(x = model, y = value)) +
    # Individual gap grey lines
    geom_line(aes(group = gap_id),
              colour = "grey75", linewidth = 0.35, alpha = 0.65) +
    geom_point(aes(colour = model), size = 1.6, alpha = 0.85) +
    # Aggregate reference black line
    geom_line(data = d_agg,
              aes(x = model, y = agg_value, group = 1),
              colour = "black", linewidth = 1.0, inherit.aes = FALSE) +
    geom_point(data = d_agg,
               aes(x = model, y = agg_value),
               colour = "black", size = 3.0, shape = 18, inherit.aes = FALSE) +
    # 4 columns × 2 rows — each ROW shares one y scale
    facet_grid(window ~ col_label, scales = "free_y") +
    scale_colour_manual(values = MODEL_COLOURS, drop = FALSE) +
    scale_x_discrete(drop = TRUE) +
    labs(x = NULL, y = metric_ylab[[metric]],
         colour = "Model",
         title  = if (add_title) title_str else NULL) +
    base_theme +
    theme(
      axis.text.x  = element_text(angle = 35, hjust = 1, size = 8),
      legend.position = "none",
      strip.text.x = element_text(size = 9,  face = "bold"),
      strip.text.y = element_text(angle = -90, size = 9, face = "bold")
    )
  
  y_vals <- c(d_indiv$value, d_agg$agg_value)
  y_vals <- y_vals[is.finite(y_vals)]
  if (length(y_vals)) {
    if (metric == "r2") {
      y_lo <- floor(min(y_vals) * 10) / 10    # round down to nearest 0.1
      y_hi <- ceiling(max(y_vals) * 10) / 10  # round up to nearest 0.1
      p <- p + coord_cartesian(ylim = c(y_lo, y_hi)) +
        scale_y_continuous(breaks = seq(y_lo, y_hi, by = 0.1))
    } else {
      p <- p + scale_y_continuous(
        breaks = seq(floor(min(y_vals)), ceiling(max(y_vals)), by = 1)
      )
    }
  }
  
  p
}

for (tv in TARGET_VARS) {
  for (gap_cat in c("S","M","L","VL")) {
    out_dir <- file.path(GRAPHS_ROOT, "site_comparison", tv, gap_cat)
    make_dir(out_dir)
    
    for (metric in c("mae","rmse","r2")) {
      p_png <- make_comparison_plot(tv, gap_cat, metric, add_title = FALSE)
      if (is.null(p_png)) next
      out_png <- file.path(out_dir,
                           paste0(toupper(metric), "_", gap_cat, ".png"))
      ggsave(out_png, p_png, width = 16, height = 7, dpi = 220, bg = "white")
      message("  Saved: ", out_png)
    }
    
    out_pdf <- file.path(out_dir,
                         paste0("comparison_", gap_cat, "_", tv, ".pdf"))
    pdf(out_pdf, width = 16, height = 7, onefile = TRUE)
    for (metric in c("mae","rmse","r2")) {
      p <- make_comparison_plot(
        tv, gap_cat, metric, add_title = TRUE,
        title_str = paste0(tv, " — ", gap_cat, " gaps — ", toupper(metric),
                           " | models ordered best\u2192worst | black \u25c6 = overall metric"))
      if (!is.null(p)) print(p)
    }
    dev.off()
    message("  Saved PDF: ", out_pdf)
  }
}


# =============================================================================
# SECTION 7 — Summary
# =============================================================================

message("\n============================================================")
message("ALL PLOTS COMPLETE")
message("============================================================")
message("Overall metrics     : ", file.path(GRAPHS_ROOT, "overall_metrics"))
message("Individual metrics  : ", file.path(GRAPHS_ROOT, "individual_metrics"))
message("Site comparison     : ", file.path(GRAPHS_ROOT, "site_comparison"))
message("Metrics CSVs        : ", csv_dir)

# ==================== end compute_metrics_and_plots.R =========================