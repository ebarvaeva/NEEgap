# =============================================================================
# METRICS AND VISUALISATION — All Models, NEE
# =============================================================================
#
# PURPOSE
#   Reads df_cv_all_predictions.rds from every model result folder, computes
#   gap-filling quality metrics (MAE, RMSE, R²) and produces:
#
#   (1) OVERALL METRICS  →  graphs/overall_metrics/NEE/
#       Line-profile plots of MAE, RMSE, R² by gap size and 30-day Grazing
#       window split.  One PNG per metric, plus a combined PDF.
#       Layout: rows = ≤30 d / >30 d window;
#               columns = JC1 Unmanaged | JC1 Managed | JC2 Unmanaged | JC2 Managed.
#
#   (2) OVERALL METRICS (NO SPLIT)  →  graphs/overall_metrics/NEE/
#       Same as (1) but pooled over ALL rows regardless of grazing window.
#       No row facet — single-row layout per metric.
#       Files: {METRIC}_overall_nosplit.png  +  appended to overall_metrics_NEE.pdf
#
# AGGREGATION METHOD
#   Overall metrics use METHOD I (pooled errors):
#     MAE(≤30d)  = sum|yhat - y| over ALL ≤30d rows across S1…Sn
#                  ─────────────────────────────────────────────────
#                  total number of ≤30d rows across S1…Sn
#   RMSE is computed analogously from pooled squared errors.
#   R² uses weighted mean of per-gap R² (weighted by n rows), since R²
#   cannot be decomposed into additive sums across gaps.
#
#   No-split metrics pool over ALL rows (le30 + gt30 combined):
#     MAE(all)  = sum|yhat - y| / n_all
#     RMSE(all) = sqrt(sum_sq_all / n_all)
#     R²(all)   = weighted.mean(r2_all, n_all)
#
# USER SETTINGS (Section 1)
#   SPLIT_DAYS controls the grazing-window boundary used in all plots.
#   Both sites (JC1, JC2) are always processed — no per-site selection needed.
#
# MODEL AVAILABILITY
#   Managed    : RF, MLP, XGBoost
#   Unmanaged  : RF, MLP, XGBoost, miniRECgap, MDS
#   MDS        : plotted as empty (placeholder) per instruction.
#
# FOLDER STRUCTURE EXPECTED
#   results/
#     {SITE}/managed/RF/df_cv_all_predictions.rds
#     {SITE}/managed/MLP/df_cv_all_predictions.rds
#     {SITE}/managed/XGBoost/df_cv_all_predictions.rds
#     {SITE}/unmanaged/RF/df_cv_all_predictions.rds
#     {SITE}/unmanaged/MLP/df_cv_all_predictions.rds
#     {SITE}/unmanaged/XGBoost/df_cv_all_predictions.rds
#     {SITE}/miniRECgap/df_cv_all_predictions.rds
#     {SITE}/MDS/df_cv_all_predictions.rds
#
# OUTPUT FOLDER STRUCTURE
#   graphs/
#     overall_metrics/
#       NEE/   — {METRIC}_overall.png         (≤30d / >30d split)
#                {METRIC}_overall_nosplit.png  (all rows, no split)
#                overall_metrics_NEE.pdf       (split pages + no-split pages)
#     metrics_csv/
#       gap_metrics_NEE_{gap_size}.csv         — per-gap, with ≤30d / >30d split
#       gap_metrics_whole_NEE_{gap_size}.csv   — per-gap, whole gap (no split)
#       overall_metrics_pooled_NEE.csv         — Method I pooled across all gaps
#                                                (one row per site × management
#                                                 × model × gap_size)
#
# =============================================================================


# =============================================================================
# SECTION 1 — User settings
# =============================================================================

# --- Grazing-window split (days) ---------------------------------------------
SPLIT_DAYS <- 30    # groups: ≤30 d  and  >30 d since last grazing event

# --- Root directories --------------------------------------------------------
RESULTS_ROOT <- here::here("results")   # where df_cv_all_predictions.rds files live
GRAPHS_ROOT  <- here::here("graphs")    # where all output plots will be saved

# --- Target variables to process ---------------------------------------------
TARGET_VARS <- c("NEE")


# =============================================================================
# SECTION 2 — Package checks and helpers
# =============================================================================

pkgs <- c("here","dplyr","tibble","tidyr","purrr","stringr","readr",
          "ggplot2","patchwork","grDevices")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tibble); library(tidyr)
  library(purrr); library(stringr); library(readr); library(ggplot2)
  library(patchwork)
})

# ---- Metric helpers (used for per-gap individual dot values) ----------------
.mae  <- function(y, yh) { ok <- is.finite(y)&is.finite(yh); if(!any(ok)) return(NA_real_); mean(abs(yh[ok]-y[ok])) }
.rmse <- function(y, yh) { ok <- is.finite(y)&is.finite(yh); if(!any(ok)) return(NA_real_); sqrt(mean((yh[ok]-y[ok])^2)) }
.r2   <- function(y, yh) { ok <- is.finite(y)&is.finite(yh); if(sum(ok)<2) return(NA_real_); cc <- suppressWarnings(cor(y[ok],yh[ok])); if(!is.finite(cc)) return(NA_real_); cc^2 }

# ---- Model catalogue --------------------------------------------------------
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

# ---- Prediction-column lookup -----------------------------------------------
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

# ---- Path builder ------------------------------------------------------------
rds_path_for_model <- function(site, model_key, management) {
  if (model_key %in% c("minirec","mds")) {
    folder <- switch(model_key, minirec="miniRECgap", mds="MDS")
    file.path(RESULTS_ROOT, site, folder, "df_cv_all_predictions.rds")
  } else {
    folder <- switch(model_key, rf="RF", mlp="MLP", xgb="XGBoost")
    file.path(RESULTS_ROOT, site, management, folder, "df_cv_all_predictions.rds")
  }
}

# ---- Safe directory creator --------------------------------------------------
make_dir <- function(...) dir.create(file.path(...), recursive=TRUE, showWarnings=FALSE)


# =============================================================================
# SECTION 3 — Load all predictions and compute metrics
# =============================================================================
# For every site × model × target × gap-size × gap label, store:
#   - Per-gap MAE/RMSE/R² for each window (used as individual dot values)
#   - Raw error SUMS and counts (used for Method I pooled aggregation)
#
# Method I pooled aggregation (overall metrics lines):
#   MAE(≤30d)  = sum_abs / n   (one division over all S1..Sn combined)
#   RMSE(≤30d) = sqrt(sum_sq / n)
#   R²(≤30d)   = weighted.mean(r2_le30, n_le30)   [cannot be pooled additively]

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
      
      flag_cols <- names(preds)[grepl(paste0("^", gap_size, "\\d+$"), names(preds))]
      if (!length(flag_cols)) next
      
      pc <- pred_col_name(target, gap_size, model_key)
      if (!pc %in% names(preds)) next
      
      truth_col <- switch(target,
                          NEE  = "NEE_orig",
                          Reco = "Reco_orig",
                          GPP  = "GPP_orig",
                          "NEE_orig")
      if (!truth_col %in% names(preds)) {
        message("  Ground-truth column '", truth_col,
                "' missing for target=", target, " in ", path,
                " — run flux_partitioning.R first. Skipping.")
        next
      }
      
      for (fc in flag_cols) {
        gap_rows <- which(preds[[fc]] %in% c(TRUE, 1))
        if (!length(gap_rows)) next
        
        y    <- as.numeric(preds[[truth_col]])[gap_rows]
        yhat <- as.numeric(preds[[pc]])[gap_rows]
        
        if (has_grazing) {
          gd     <- as.numeric(preds$Grazing_days_since)[gap_rows]
          sel_le <- is.finite(gd) & gd <= SPLIT_DAYS
          sel_gt <- is.finite(gd) & gd >  SPLIT_DAYS
        } else {
          sel_le <- rep(TRUE,  length(gap_rows))
          sel_gt <- rep(FALSE, length(gap_rows))
        }
        
        # Valid indices per window (finite on both y and yhat)
        ok_all <- is.finite(y) & is.finite(yhat)          # all gap rows, no split
        ok_le  <- sel_le & ok_all
        ok_gt  <- sel_gt & ok_all
        
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
          
          # Per-gap metrics — used as individual dot values in comparison plots
          mae_all    = .mae(y, yhat),
          rmse_all   = .rmse(y, yhat),
          r2_all     = .r2(y, yhat),
          mae_le30   = .mae(y[sel_le], yhat[sel_le]),
          rmse_le30  = .rmse(y[sel_le], yhat[sel_le]),
          r2_le30    = .r2(y[sel_le], yhat[sel_le]),
          mae_gt30   = .mae(y[sel_gt], yhat[sel_gt]),
          rmse_gt30  = .rmse(y[sel_gt], yhat[sel_gt]),
          r2_gt30    = .r2(y[sel_gt], yhat[sel_gt]),
          
          # Raw error sums — used for Method I pooled aggregation
          # MAE_pooled  = sum_abs / n   (one division over all S1..Sn combined)
          # RMSE_pooled = sqrt(sum_sq / n)
          n_all        = sum(ok_all),
          sum_abs_all  = sum(abs(yhat[ok_all] - y[ok_all])),
          sum_sq_all   = sum((yhat[ok_all]    - y[ok_all])^2),
          n_le30       = sum(ok_le),
          n_gt30       = sum(ok_gt),
          sum_abs_le30 = sum(abs(yhat[ok_le] - y[ok_le])),
          sum_sq_le30  = sum((yhat[ok_le]    - y[ok_le])^2),
          sum_abs_gt30 = sum(abs(yhat[ok_gt] - y[ok_gt])),
          sum_sq_gt30  = sum((yhat[ok_gt]    - y[ok_gt])^2)
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
    model      = factor(model, levels = c("RF","MLP","XGBoost","miniRECgap","MDS")),
    gap_size   = factor(gap_size, levels = c("S","M","L","VL")),
    management = factor(management, levels = c("unmanaged","managed"))
  )

message("Gap-level metrics computed: ", nrow(all_gap_metrics), " rows.")

# ---- Save CSVs — per-gap metrics with ≤30d / >30d split --------------------
# One file per target × gap size.  Each row = one gap Si.
# Columns include: mae_le30, rmse_le30, r2_le30, mae_gt30, rmse_gt30, r2_gt30,
#                  n_le30, n_gt30, and the raw sums used for Method I pooling.
csv_dir <- file.path(GRAPHS_ROOT, "metrics_csv")
make_dir(csv_dir)
for (tv in TARGET_VARS) {
  for (gs in c("S","M","L","VL")) {
    d <- filter(all_gap_metrics, target == tv, gap_size == gs)
    write_csv(d, file.path(csv_dir, paste0("gap_metrics_", tv, "_", gs, ".csv")))
  }
}
message("Per-gap split CSVs saved to: ", csv_dir)

# ---- Save CSVs — per-gap metrics, no 30d split (whole gap Si) ---------------
# One file per target × gap size.  Each row = one gap Si.
# Columns: mae, rmse, r2 — computed on ALL rows of that gap regardless of
# Grazing_days_since.  Simpler summary for analyses that do not need the split.
for (tv in TARGET_VARS) {
  for (gs in c("S","M","L","VL")) {
    d <- filter(all_gap_metrics, target == tv, gap_size == gs) %>%
      select(site, management, model, model_key, target, gap_size, gap_id, gap_num,
             mae  = mae_all,
             rmse = rmse_all,
             r2   = r2_all)
    write_csv(d, file.path(csv_dir, paste0("gap_metrics_whole_", tv, "_", gs, ".csv")))
  }
}
message("Whole-gap CSVs saved to: ", csv_dir)

# ---- Drop MDS from VL (MDS cannot fill 30-day gaps) ------------------------
all_gap_metrics <- all_gap_metrics %>%
  filter(!(model == "MDS" & gap_size %in% c("L", "VL")))

# ---- Overall metrics — METHOD I pooled aggregation -------------------------
# MAE  = sum of all absolute errors across S1…Sn  /  total number of points
# RMSE = sqrt(sum of squared errors / total points)
# R²   = weighted mean of per-gap R² (weighted by n, since R² is not additive)
overall_metrics <- all_gap_metrics %>%
  group_by(site, management, model, model_key, target, gap_size) %>%
  summarise(
    mae_all   = sum(sum_abs_all, na.rm=TRUE) / max(sum(n_all, na.rm=TRUE), 1),
    rmse_all  = sqrt(sum(sum_sq_all, na.rm=TRUE) / max(sum(n_all, na.rm=TRUE), 1)),
    r2_all    = weighted.mean(r2_all, w = pmax(n_all, 0), na.rm=TRUE),
    mae_le30  = sum(sum_abs_le30, na.rm=TRUE) / max(sum(n_le30, na.rm=TRUE), 1),
    rmse_le30 = sqrt(sum(sum_sq_le30, na.rm=TRUE) / max(sum(n_le30, na.rm=TRUE), 1)),
    r2_le30   = weighted.mean(r2_le30, w = pmax(n_le30, 0), na.rm=TRUE),
    mae_gt30  = sum(sum_abs_gt30, na.rm=TRUE) / max(sum(n_gt30, na.rm=TRUE), 1),
    rmse_gt30 = sqrt(sum(sum_sq_gt30, na.rm=TRUE) / max(sum(n_gt30, na.rm=TRUE), 1)),
    r2_gt30   = weighted.mean(r2_gt30, w = pmax(n_gt30, 0), na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .)))

# ---- Save CSVs — pooled Method I metrics (one row per model × gap size) -----
# These are the exact values shown as lines in the overall_metrics plots.
# Columns without split: mae_all, rmse_all, r2_all  — pooled over ALL gap rows
# Columns with split   : mae_le30, rmse_le30, r2_le30, mae_gt30, rmse_gt30, r2_gt30
# MAE/RMSE use Method I pooling (sum errors / total n).
# R² uses weighted mean of per-gap R² (weighted by n).
for (tv in TARGET_VARS) {
  write_csv(
    filter(overall_metrics, target == tv),
    file.path(csv_dir, paste0("overall_metrics_pooled_", tv, ".csv"))
  )
}
message("Pooled overall metrics CSVs saved to: ", csv_dir)


# =============================================================================
# SECTION 4 — Plot theme and shared helpers
# =============================================================================

base_theme <- theme_bw(base_size = 15) +
  theme(
    panel.grid.minor    = element_blank(),
    strip.text          = element_text(size = 14),
    legend.position     = "bottom",
    legend.title        = element_text(size = 14),
    legend.text         = element_text(size = 14),
    plot.background     = element_rect(fill = "white", colour = NA),
    panel.spacing       = unit(1, "lines")
  )

MODEL_COLOURS <- c(
  RF          = "#E41A1C",
  MLP         = "#377EB8",
  XGBoost     = "#4DAF4A",
  miniRECgap  = "#FF7F00",
  MDS         = "#999999"
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

COL_LABEL_LEVELS <- c("JC1 Unmanaged","JC1 Managed","JC2 Unmanaged","JC2 Managed")


# =============================================================================
# SECTION 5 — Overall metrics plots  (JC1 + JC2, managed + unmanaged)
# =============================================================================

message("\nBuilding OVERALL METRICS plots (4-column: JC1/JC2 × managed/unmanaged) ...")

# ---- 5a: Split plot (≤30d / >30d grazing window rows) ----------------------
# Layout: rows = window (le30 / gt30), columns = site × management
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
      strip.text.y    = element_text(angle = -90, size = 14)
    )
  
  y_vals <- d$value[is.finite(d$value)]
  if (length(y_vals)) {
    if (metric == "r2") {
      y_lo <- floor(min(y_vals) * 10) / 10
      y_hi <- ceiling(max(y_vals) * 10) / 10
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

# ---- 5b: No-split plot (all rows pooled, no grazing-window facet row) -------
# Layout: single row, columns = site × management
# Uses _all columns from overall_metrics (MAE/RMSE pooled over le30 + gt30).
make_overall_plot_nosplit <- function(tv, metric, add_title = FALSE, title_str = "") {
  
  d <- overall_metrics %>%
    filter(target == tv) %>%
    mutate(
      col_label = factor(paste(site, str_to_title(management)),
                         levels = COL_LABEL_LEVELS)
    ) %>%
    select(col_label, model, gap_size,
           value = paste0(metric, "_all")) %>%
    filter(!is.na(value))
  
  if (!nrow(d)) return(NULL)
  
  p <- ggplot(d, aes(x = gap_size, y = value,
                     colour = model, group = model)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    facet_grid(. ~ col_label, scales = "free_y") +   # single row, no window facet
    scale_colour_manual(values = MODEL_COLOURS, drop = FALSE) +
    labs(x = "Gap size", y = metric_ylab[[metric]],
         colour = "Model",
         title  = if (add_title) title_str else NULL) +
    base_theme +
    theme(legend.position = "bottom")
  
  y_vals <- d$value[is.finite(d$value)]
  if (length(y_vals)) {
    if (metric == "r2") {
      y_lo <- floor(min(y_vals) * 10) / 10
      y_hi <- ceiling(max(y_vals) * 10) / 10
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

# ---- Generate and save all plots --------------------------------------------
for (tv in TARGET_VARS) {
  out_dir <- file.path(GRAPHS_ROOT, "overall_metrics", tv)
  make_dir(out_dir)
  
  pdf_plots <- list()
  
  for (metric in c("mae","rmse","r2")) {
    
    # -- Split plot (≤30d / >30d) ---
    p_split <- make_overall_plot(tv, metric, add_title = FALSE)
    if (!is.null(p_split)) {
      out_png <- file.path(out_dir, paste0(toupper(metric), "_overall.png"))
      ggsave(out_png, p_split, width = 14, height = 7, dpi = 220, bg = "white")
      message("  Saved: ", out_png)
      pdf_plots[[paste0(metric, "_split")]] <- make_overall_plot(
        tv, metric, add_title = TRUE,
        title_str = paste0(tv, " — ", toupper(metric),
                           " by gap size | JC1 & JC2, managed & unmanaged",
                           " | pooled MAE/RMSE across all gaps in window"))
    }
    
    # -- No-split plot (all rows) ---
    p_nosplit <- make_overall_plot_nosplit(tv, metric, add_title = FALSE)
    if (!is.null(p_nosplit)) {
      out_png_ns <- file.path(out_dir, paste0(toupper(metric), "_overall_nosplit.png"))
      ggsave(out_png_ns, p_nosplit, width = 14, height = 4, dpi = 220, bg = "white")
      message("  Saved: ", out_png_ns)
      pdf_plots[[paste0(metric, "_nosplit")]] <- make_overall_plot_nosplit(
        tv, metric, add_title = TRUE,
        title_str = paste0(tv, " — ", toupper(metric),
                           " by gap size | JC1 & JC2, managed & unmanaged",
                           " | pooled over ALL rows (no grazing-window split)"))
    }
  }
  
  # PDF: split pages first (mae, rmse, r2), then no-split pages
  out_pdf <- file.path(out_dir, paste0("overall_metrics_", tv, ".pdf"))
  pdf(out_pdf, width = 14, height = 7, onefile = TRUE)
  for (key in c("mae_split","rmse_split","r2_split",
                "mae_nosplit","rmse_nosplit","r2_nosplit")) {
    p <- pdf_plots[[key]]
    if (!is.null(p)) print(p)
  }
  dev.off()
  message("  Saved PDF: ", out_pdf)
}


# =============================================================================
# SECTION 6 — Summary
# =============================================================================

message("\n============================================================")
message("ALL PLOTS COMPLETE")
message("============================================================")
message("Overall metrics  : ", file.path(GRAPHS_ROOT, "overall_metrics"))
message("  Split PNGs     : {METRIC}_overall.png        (≤30d / >30d rows)")
message("  No-split PNGs  : {METRIC}_overall_nosplit.png (all rows pooled)")
message("  Combined PDF   : overall_metrics_{tv}.pdf")
message("Metrics CSVs     : ", csv_dir)

# ==================== end compute_metrics_and_plots.R =========================