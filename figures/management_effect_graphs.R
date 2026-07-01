# management_effect_graphs.R — Management variable contribution: CSVs + % MAE-reduction heatmaps
#
# Quantifies how much each management predictor helps NEE gap-filling by adding it
# one at a time on top of the meteorology-only baseline (BASE) and measuring the
# out-of-gap MAE change:
#   ΔMAE_m      = MAE(BASE+m) − MAE(BASE)          (negative = improvement)
#   Reduction_m = −100 × ΔMAE_m / MAE(BASE)  (%)   (positive = improvement)
# Metrics use Method I pooling (identical to compute_metrics_and_plots.R): one
# division over all half-hours in a gap-size class, and R² is an n-weighted mean of
# per-gap R². Only RF/MLP/XGBoost are used (MDS and miniRECgap have fixed predictor
# structures and cannot take extra covariates).
#
# Inputs:
#   results/{SITE}/management_effect/{MODEL}/{MGMT_VAR}/df_cv_all_predictions.rds
#   graphs/metrics_csv/overall_metrics_pooled_NEE.csv  (BASE, from
#       compute_metrics_and_plots.R; must carry site + management columns —
#       re-run that script if they are absent)
#
# Management variables (m), with the heatmap y-axis symbol:
#   dsg   Grazing_days_since     — days since last grazing
#   dsf   Fertiliser_days_since  — days since last fertilisation
#   h[s]  grass_height           — sward height
#   B     grass_biomass          — above-ground dry-matter biomass
#   PI    phytomass_index        — flux-derived canopy-state index
#
# Outputs:
#   graphs/metrics_csv/management_effect/{SITE}/{MGMT_VAR}/overall_metrics_pooled_NEE.csv
#   graphs/management_effect/heatmaps/NEE/pct_MAE_{S,M,L,VL}.png
#       (green = improvement over BASE, red = degradation; pooled over all gap rows)
#
# To adapt: change SITES / MODELS / MGMT_VARS (and MGMT_LABELS), or GAP_SIZES.


# Repository roots and the heatmap output root
RESULTS_ROOT <- here::here("results")
GRAPHS_ROOT  <- here::here("graphs")
OUT_ROOT     <- file.path(GRAPHS_ROOT, "management_effect")

# Sites, models, target flux, gap-size classes
SITES       <- c("JC1", "JC2")
MODELS      <- c("RF", "MLP", "XGBoost")
TARGET_VARS <- c("NEE")
GAP_SIZES   <- c("S", "M", "L", "VL")
SPLIT_DAYS  <- 30   # kept for CSV completeness; not used in plots

# The five management predictors tested, one per ablation run
MGMT_VARS <- c(
  "Grazing_days_since",
  "Fertiliser_days_since",
  "grass_height",
  "grass_biomass",
  "phytomass_index"
)

# Heatmap y-axis symbols; subscript strings (e.g. "h[s]") are parsed as math text
MGMT_LABELS <- c(
  Grazing_days_since    = "dsg",
  Fertiliser_days_since = "dsf",
  grass_height          = "h[s]",
  grass_biomass         = "B",
  phytomass_index       = "PI"
)

# Prediction-column suffix per model inside each df_cv_all_predictions.rds
PRED_SUFFIX <- c(RF = "rf_predicted", MLP = "mlp_predicted", XGBoost = "xgb_predicted")


# Libraries
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(readr); library(ggplot2); library(tibble)
  library(forcats); library(lubridate)
})

# Create a (possibly nested) output directory if it does not exist
make_dir <- function(...) dir.create(file.path(...), recursive = TRUE, showWarnings = FALSE)

# Method I metric helpers — identical to compute_metrics_and_plots.R
.mae  <- function(y, yh) { ok <- is.finite(y) & is.finite(yh); if (!any(ok)) return(NA_real_); mean(abs(yh[ok] - y[ok])) }
.rmse <- function(y, yh) { ok <- is.finite(y) & is.finite(yh); if (!any(ok)) return(NA_real_); sqrt(mean((yh[ok] - y[ok])^2)) }
.r2   <- function(y, yh) { ok <- is.finite(y) & is.finite(yh); if (sum(ok) < 2) return(NA_real_); cc <- suppressWarnings(cor(y[ok], yh[ok])); if (!is.finite(cc)) return(NA_real_); cc^2 }

# Shared heatmap theme
base_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(size = 9),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 8),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.spacing    = unit(0.8, "lines")
  )


# Read one management-effect RDS and return per-gap raw sums + per-gap R² (Method I inputs)
compute_mgmt_gap_metrics <- function(site, model, mgmt_var) {
  rds_path <- file.path(
    RESULTS_ROOT, site, "management_effect", model, mgmt_var,
    "df_cv_all_predictions.rds"
  )
  if (!file.exists(rds_path)) { return(tibble()) }

  preds <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(preds)) { return(tibble()) }

  suffix   <- PRED_SUFFIX[[model]]
  has_graz <- "Grazing_days_since" %in% names(preds)

  rows <- list()

  for (target in TARGET_VARS) {
    truth_col <- switch(target, NEE = "NEE_orig", Reco = "Reco_orig", GPP = "GPP_orig")
    if (!truth_col %in% names(preds)) next

    for (gs in GAP_SIZES) {
      pc        <- paste0(target, "_", gs, "_", suffix)
      flag_cols <- names(preds)[grepl(paste0("^", gs, "\\d+$"), names(preds))]
      if (!pc %in% names(preds) || !length(flag_cols)) next

      # One row per artificial gap of this size
      for (fc in flag_cols) {
        gap_rows <- which(preds[[fc]] %in% c(TRUE, 1L))
        if (!length(gap_rows)) next

        y    <- as.numeric(preds[[truth_col]])[gap_rows]
        yhat <- as.numeric(preds[[pc]])[gap_rows]
        sub  <- preds[gap_rows, ]

        # Window splits — computed for CSV completeness, not used in plots
        if (has_graz) {
          gd     <- as.numeric(sub$Grazing_days_since)
          sel_le <- is.finite(gd) & gd <= SPLIT_DAYS
          sel_gt <- is.finite(gd) & gd >  SPLIT_DAYS
        } else {
          sel_le <- rep(TRUE,  length(gap_rows))
          sel_gt <- rep(FALSE, length(gap_rows))
        }

        ok_all <- is.finite(y) & is.finite(yhat)
        ok_le  <- sel_le & ok_all
        ok_gt  <- sel_gt & ok_all

        gap_num <- as.integer(str_extract(fc, "\\d+$"))

        rows[[length(rows) + 1]] <- tibble(
          site       = site,
          management = paste0("BASE+", mgmt_var),
          model      = model,
          model_key  = tolower(str_extract(model, "[A-Za-z]+")),
          mgmt_var   = mgmt_var,
          target     = target,
          gap_size   = gs,
          gap_id     = paste0(gs, gap_num),
          gap_num    = gap_num,
          # Per-gap R² for weighted-mean pooling
          r2_all     = .r2(y,          yhat),
          r2_le30    = .r2(y[sel_le],  yhat[sel_le]),
          r2_gt30    = .r2(y[sel_gt],  yhat[sel_gt]),
          # Raw sums for Method I MAE/RMSE pooling
          n_all           = sum(ok_all),
          sum_abs_all     = sum(abs(yhat[ok_all] - y[ok_all])),
          sum_sq_all      = sum((yhat[ok_all]    - y[ok_all])^2),
          n_le30          = sum(ok_le),
          sum_abs_le30    = sum(abs(yhat[ok_le]  - y[ok_le])),
          sum_sq_le30     = sum((yhat[ok_le]      - y[ok_le])^2),
          n_gt30          = sum(ok_gt),
          sum_abs_gt30    = sum(abs(yhat[ok_gt]  - y[ok_gt])),
          sum_sq_gt30     = sum((yhat[ok_gt]      - y[ok_gt])^2)
        )
      }
    }
  }
  if (!length(rows)) return(tibble())
  bind_rows(rows)
}

# Gap-level metrics for every site × model × management-variable run
all_mgmt_gap <- map_dfr(SITES, function(site) {
  map_dfr(MODELS, function(model) {
    map_dfr(MGMT_VARS, function(mv) {
      compute_mgmt_gap_metrics(site, model, mv)
    })
  })
}) %>%
  mutate(gap_size = factor(gap_size, levels = GAP_SIZES))

# Method I pooled aggregation over all gaps within each size class
pooled_mgmt <- all_mgmt_gap %>%
  group_by(site, management, model, model_key, mgmt_var, target, gap_size) %>%
  summarise(
    mae_all   = sum(sum_abs_all,  na.rm = TRUE) / max(sum(n_all,  na.rm = TRUE), 1),
    rmse_all  = sqrt(sum(sum_sq_all,  na.rm = TRUE) / max(sum(n_all,  na.rm = TRUE), 1)),
    r2_all    = weighted.mean(r2_all,  w = pmax(n_all,  0), na.rm = TRUE),
    mae_le30  = sum(sum_abs_le30, na.rm = TRUE) / max(sum(n_le30, na.rm = TRUE), 1),
    rmse_le30 = sqrt(sum(sum_sq_le30, na.rm = TRUE) / max(sum(n_le30, na.rm = TRUE), 1)),
    r2_le30   = weighted.mean(r2_le30, w = pmax(n_le30, 0), na.rm = TRUE),
    mae_gt30  = sum(sum_abs_gt30, na.rm = TRUE) / max(sum(n_gt30, na.rm = TRUE), 1),
    rmse_gt30 = sqrt(sum(sum_sq_gt30, na.rm = TRUE) / max(sum(n_gt30, na.rm = TRUE), 1)),
    r2_gt30   = weighted.mean(r2_gt30, w = pmax(n_gt30, 0), na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA_real_, .)))

# Write one pooled CSV per site × management variable
for (site in SITES) {
  for (mv in MGMT_VARS) {
    csv_out_dir <- file.path(GRAPHS_ROOT, "metrics_csv", "management_effect", site, mv)
    make_dir(csv_out_dir)
    for (tv in TARGET_VARS) {
      d <- pooled_mgmt %>%
        filter(site == .env$site, mgmt_var == mv, target == tv) %>%
        select(site, management, model, model_key, target, gap_size,
               mae_all, rmse_all, r2_all,
               mae_le30, rmse_le30, r2_le30,
               mae_gt30, rmse_gt30, r2_gt30)
      write_csv(d, file.path(csv_out_dir, paste0("overall_metrics_pooled_", tv, ".csv")))
    }
  }
}


# Load the BASE (unmanaged) pooled metrics used as the delta reference
base_csv_dir <- file.path(GRAPHS_ROOT, "metrics_csv")

base_pooled <- map_dfr(TARGET_VARS, function(tv) {
  p <- file.path(base_csv_dir, paste0("overall_metrics_pooled_", tv, ".csv"))
  if (!file.exists(p)) { return(tibble()) }
  d <- read_csv(p, show_col_types = FALSE)
  # Fail early if the BASE CSV lacks the columns needed to match runs
  required <- c("site", "management", "model", "target", "gap_size", "mae_all")
  missing  <- setdiff(required, names(d))
  if (length(missing))
    stop("overall_metrics_pooled_", tv, ".csv is missing columns: ",
         paste(missing, collapse = ", "),
         "\nRe-run compute_metrics_and_plots.R to regenerate with site + management columns.")
  d
}) %>%
  filter(management == "unmanaged", model %in% MODELS) %>%
  select(site, model, target, gap_size,
         base_mae_all  = mae_all,
         base_rmse_all = rmse_all,
         base_r2_all   = r2_all)


# Join each BASE+m run to its BASE and compute deltas + % MAE reduction
delta_pooled <- pooled_mgmt %>%
  left_join(base_pooled, by = c("site", "model", "target", "gap_size")) %>%
  mutate(
    delta_mae_all     = mae_all  - base_mae_all,
    delta_rmse_all    = rmse_all - base_rmse_all,
    delta_r2_all      = r2_all   - base_r2_all,
    pct_mae_reduction = -100 * delta_mae_all / base_mae_all,
    # mgmt_label uses symbol strings parsed as math text in heatmap y-axis
    mgmt_label        = factor(MGMT_LABELS[mgmt_var], levels = unname(MGMT_LABELS)),
    gap_size          = factor(gap_size, levels = GAP_SIZES)
  )


# Build one % MAE-reduction heatmap (variable × model, faceted by site) for a gap size
make_heatmap <- function(tv, gs) {
  d <- delta_pooled %>%
    filter(target == tv, gap_size == gs, !is.na(pct_mae_reduction)) %>%
    mutate(model = factor(model, levels = MODELS))
  if (!nrow(d)) return(NULL)

  # Symmetric colour limit rounded up to the next multiple of 5 %
  lim <- ceiling(max(abs(d$pct_mae_reduction), na.rm = TRUE) / 5) * 5

  ggplot(d, aes(x = model, y = mgmt_label, fill = pct_mae_reduction)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%+.1f%%", pct_mae_reduction)), size = 3) +
    facet_wrap(~site, nrow = 1) +
    scale_fill_gradient2(
      low = "#D73027", mid = "white", high = "#1A9850",
      midpoint = 0, limits = c(-lim, lim),
      name = "% MAE\nreduction"
    ) +
    # Parse y-axis labels as math expressions so subscripts render correctly
    scale_y_discrete(labels = function(x) parse(text = x)) +
    labs(
      x        = NULL,
      y        = NULL,
      subtitle = paste0(tv, " \u2014 ", gs, " gaps")
    ) +
    base_theme +
    theme(legend.position = "right")
}

# One heatmap PNG per gap size
for (tv in TARGET_VARS) {
  out_dir <- file.path(OUT_ROOT, "heatmaps", tv)
  make_dir(out_dir)
  for (gs in GAP_SIZES) {
    p <- make_heatmap(tv, gs)
    if (is.null(p)) next
    ggsave(
      file.path(out_dir, paste0("pct_MAE_", gs, ".png")),
      p, width = 9, height = 5, dpi = 220, bg = "white"
    )
  }
}
