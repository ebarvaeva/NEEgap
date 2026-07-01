# compute_metrics_and_plots.R — Gap-filling metrics + overall-metric plots (run first)
#
# Reads df_cv_all_predictions.rds from every model result folder, computes per-gap and pooled MAE / RMSE / R2
# for NEE, writes the CSVs that the other metrics scripts consume, and draws the
# overall-metric line plots (gap size on x, model as colour).
#
# Pooling (Method I): MAE = sum|yhat-y| / N and RMSE = sqrt(sum(yhat-y)^2 / N)
# over all half-hours in a gap-size category; R2 = n-weighted mean of per-gap R2
# (R2 is not additive across gaps). "le30" / "gt30" restrict the sum to rows
# <=30 d / >30 d since the last grazing event; "all" pools both.
#
# Inputs (per site, under results/{SITE}/):
#   managed/{RF,MLP,XGBoost}/df_cv_all_predictions.rds
#   unmanaged/{RF,MLP,XGBoost}/df_cv_all_predictions.rds
#   miniRECgap/df_cv_all_predictions.rds
#   MDS/df_cv_all_predictions.rds        (MDS excluded from L and VL gaps)
#   each needs: NEE_orig, gap-flag columns S1.., M1.., L1.., VL1.., prediction
#   columns NEE_{S|M|L|VL}_{rf|mlp|xgb|minirec|mds}_predicted, and (optional)
#   Grazing_days_since for the window split.
#
# Outputs:
#   graphs/metrics_csv/gap_metrics_NEE_{S|M|L|VL}.csv          (per-gap, split)
#   graphs/metrics_csv/gap_metrics_whole_NEE_{S|M|L|VL}.csv    (per-gap, no split)
#   graphs/metrics_csv/overall_metrics_pooled_NEE.csv          (pooled, all models)
#   graphs/overall_metrics/NEE/{MAE|RMSE|R2}_overall.png       (le30 / gt30 rows)
#   graphs/overall_metrics/NEE/{MAE|RMSE|R2}_overall_nosplit.png
#   graphs/overall_metrics/NEE/overall_metrics_NEE.pdf
#
# To adapt: add sites to SITES_ALL; set SPLIT_DAYS to change the grazing window.

# Grazing-window split (days): rows grouped as <=SPLIT_DAYS and >SPLIT_DAYS.
SPLIT_DAYS <- 30

# Root directories.
RESULTS_ROOT <- here::here("results")   # where df_cv_all_predictions.rds files live
GRAPHS_ROOT  <- here::here("graphs")    # where all output plots will be saved

TARGET_VARS <- c("NEE")

# Fail early if any required package is missing.
pkgs <- c("here","dplyr","tibble","tidyr","purrr","stringr","readr",
          "ggplot2","patchwork","grDevices")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install missing packages: ", paste(miss, collapse = ", "))

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tibble); library(tidyr)
  library(purrr); library(stringr); library(readr); library(ggplot2)
  library(patchwork)
})

# Per-gap metric helpers (finite pairs only); used for individual-gap dot values.
.mae  <- function(y, yh) { ok <- is.finite(y)&is.finite(yh); if(!any(ok)) return(NA_real_); mean(abs(yh[ok]-y[ok])) }
.rmse <- function(y, yh) { ok <- is.finite(y)&is.finite(yh); if(!any(ok)) return(NA_real_); sqrt(mean((yh[ok]-y[ok])^2)) }
.r2   <- function(y, yh) { ok <- is.finite(y)&is.finite(yh); if(sum(ok)<2) return(NA_real_); cc <- suppressWarnings(cor(y[ok],yh[ok])); if(!is.finite(cc)) return(NA_real_); cc^2 }

# Model catalogue: which model/management combinations to look for.
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

# pred_col_name(): prediction column for a target x gap-size x model.
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

# rds_path_for_model(): result-file path (miniRECgap/MDS sit outside managed/unmanaged).
rds_path_for_model <- function(site, model_key, management) {
  if (model_key %in% c("minirec","mds")) {
    folder <- switch(model_key, minirec="miniRECgap", mds="MDS")
    file.path(RESULTS_ROOT, site, folder, "df_cv_all_predictions.rds")
  } else {
    folder <- switch(model_key, rf="RF", mlp="MLP", xgb="XGBoost")
    file.path(RESULTS_ROOT, site, management, folder, "df_cv_all_predictions.rds")
  }
}

make_dir <- function(...) dir.create(file.path(...), recursive=TRUE, showWarnings=FALSE)

# compute_metrics_one_model(): for one site x model, return one row per gap Si with
# per-gap MAE/RMSE/R2 (dot values) and raw error sums + counts (Method I pooling).
compute_metrics_one_model <- function(site, model_key, management) {
  
  path <- rds_path_for_model(site, model_key, management)
  if (!file.exists(path)) {
    return(tibble())
  }
  
  preds <- tryCatch(readRDS(path), error = function(e) {
    return(NULL)
  })
  if (is.null(preds)) return(tibble())
  
  if (!"NEE_orig" %in% names(preds)) {
    return(tibble())
  }
  
  has_grazing <- "Grazing_days_since" %in% names(preds)
  
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
      
      # Ground-truth column for the target (NEE here; Reco/GPP handled if added).
      truth_col <- switch(target,
                          NEE  = "NEE_orig",
                          Reco = "Reco_orig",
                          GPP  = "GPP_orig",
                          "NEE_orig")
      if (!truth_col %in% names(preds)) {
        next
      }
      
      for (fc in flag_cols) {
        gap_rows <- which(preds[[fc]] %in% c(TRUE, 1))
        if (!length(gap_rows)) next
        
        y    <- as.numeric(preds[[truth_col]])[gap_rows]
        yhat <- as.numeric(preds[[pc]])[gap_rows]
        
        # Window masks: <=SPLIT_DAYS vs >SPLIT_DAYS since last grazing.
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

# Run over every site x model combination and factor the key columns.
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

# Per-gap split CSVs (one row per gap Si, with le30 / gt30 columns and raw sums).
csv_dir <- file.path(GRAPHS_ROOT, "metrics_csv")
make_dir(csv_dir)
for (tv in TARGET_VARS) {
  for (gs in c("S","M","L","VL")) {
    d <- filter(all_gap_metrics, target == tv, gap_size == gs)
    write_csv(d, file.path(csv_dir, paste0("gap_metrics_", tv, "_", gs, ".csv")))
  }
}

# Per-gap whole-gap CSVs (same rows, MAE/RMSE/R2 over all rows, no window split).
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

# MDS cannot fill L/VL gaps — drop those rows.
all_gap_metrics <- all_gap_metrics %>%
  filter(!(model == "MDS" & gap_size %in% c("L", "VL")))

# Method I pooled aggregation (these values are the lines in the overall plots).
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

# Pooled Method I CSVs (one row per site x management x model x gap size).
for (tv in TARGET_VARS) {
  write_csv(
    filter(overall_metrics, target == tv),
    file.path(csv_dir, paste0("overall_metrics_pooled_", tv, ".csv"))
  )
}

# Shared plot theme and label lookups.
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

# make_overall_plot(): split plot — rows = le30 / gt30 window, columns = site x management.
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

# make_overall_plot_nosplit(): single row, columns = site x management, pooled over all rows.
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

# For each metric: save split PNG, no-split PNG, and collect titled versions for the PDF.
for (tv in TARGET_VARS) {
  out_dir <- file.path(GRAPHS_ROOT, "overall_metrics", tv)
  make_dir(out_dir)
  
  pdf_plots <- list()
  
  for (metric in c("mae","rmse","r2")) {
    
    p_split <- make_overall_plot(tv, metric, add_title = FALSE)
    if (!is.null(p_split)) {
      out_png <- file.path(out_dir, paste0(toupper(metric), "_overall.png"))
      ggsave(out_png, p_split, width = 14, height = 7, dpi = 220, bg = "white")
      pdf_plots[[paste0(metric, "_split")]] <- make_overall_plot(
        tv, metric, add_title = TRUE,
        title_str = paste0(tv, " — ", toupper(metric),
                           " by gap size | JC1 & JC2, managed & unmanaged",
                           " | pooled MAE/RMSE across all gaps in window"))
    }
    
    p_nosplit <- make_overall_plot_nosplit(tv, metric, add_title = FALSE)
    if (!is.null(p_nosplit)) {
      out_png_ns <- file.path(out_dir, paste0(toupper(metric), "_overall_nosplit.png"))
      ggsave(out_png_ns, p_nosplit, width = 14, height = 4, dpi = 220, bg = "white")
      pdf_plots[[paste0(metric, "_nosplit")]] <- make_overall_plot_nosplit(
        tv, metric, add_title = TRUE,
        title_str = paste0(tv, " — ", toupper(metric),
                           " by gap size | JC1 & JC2, managed & unmanaged",
                           " | pooled over ALL rows (no grazing-window split)"))
    }
  }
  
  # PDF: split pages first (mae, rmse, r2), then no-split pages.
  out_pdf <- file.path(out_dir, paste0("overall_metrics_", tv, ".pdf"))
  pdf(out_pdf, width = 14, height = 7, onefile = TRUE)
  for (key in c("mae_split","rmse_split","r2_split",
                "mae_nosplit","rmse_nosplit","r2_nosplit")) {
    p <- pdf_plots[[key]]
    if (!is.null(p)) print(p)
  }
  dev.off()
}
