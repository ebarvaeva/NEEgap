# =============================================================================
# MANAGEMENT EFFECT ANALYSIS
# Feature-Addition Ablation Study — CSV Export & Heatmaps
# =============================================================================
#
# PURPOSE
#   Quantifies the marginal contribution of individual management predictors
#   to NEE gap-filling performance by comparing each BASE+m run against the
#   meteorology-only baseline (BASE).
#
#   For each management variable m added one at a time on top of BASE:
#     ΔMAE_m        = MAE(BASE+m) − MAE(BASE)          [negative = improvement]
#     Reduction_m   = −100 × ΔMAE_m / MAE(BASE)  (%)  [positive = improvement]
#
# INPUTS
#   BASE predictions (meteorological + temporal predictors only):
#     results/{SITE}/unmanaged/{MODEL}/df_cv_all_predictions.rds
#
#   Management-effect predictions (BASE + one variable at a time):
#     results/{SITE}/management_effect/{MODEL}/{MGMT_VAR}/df_cv_all_predictions.rds
#
#   BASE pooled metrics CSV (written by compute_metrics_and_plots.R):
#     graphs/metrics_csv/overall_metrics_pooled_NEE.csv
#     Requires columns: site, management, model, target, gap_size, mae_all
#     → used directly for heatmap delta computation
#     → if site/management columns are missing, re-run compute_metrics_and_plots.R
#
# OUTPUTS
#   1. Management-effect pooled CSVs  (Method I — same pooling as BASE CSVs)
#        graphs/metrics_csv/management_effect/{SITE}/{MGMT_VAR}/
#          overall_metrics_pooled_NEE.csv
#      Columns: site, management (= "BASE+{var}"), model, model_key, target,
#               gap_size, mae_all, rmse_all, r2_all,
#               mae_le30, rmse_le30, r2_le30,   <- kept for completeness
#               mae_gt30, rmse_gt30, r2_gt30    <- kept for completeness
#
#   2. Heatmaps  —  % MAE reduction per variable x model, faceted by site
#        graphs/management_effect/heatmaps/NEE/pct_MAE_{GAP_SIZE}.png
#      One PNG per gap size (S / M / L / VL).
#      Metrics pooled over ALL gap rows (no grazing-window split).
#      y-axis labels use table symbols (dsg, dsf, h_s, B, N, PI).
#      Green = improvement over BASE; red = degradation.
#
# MODELS INCLUDED
#   RF, MLP, XGBoost  (MDS and miniRECgap excluded — their predictor
#   structures are fixed and cannot be extended with arbitrary covariates)
#
# MANAGEMENT VARIABLES (m) — symbol as in Table tab:variables
#   dsg  (Grazing_days_since)    — days since last grazing event
#   dsf  (Fertiliser_days_since) — days since last fertilisation event
#   N    (N)                     — nitrogen application amount
#   h_s  (grass_height)          — sward height
#   B    (grass_biomass)         — above-ground dry-matter biomass
#   PI   (phytomass_index)       — composite canopy state index
#
# AGGREGATION METHOD  (Method I — consistent with compute_metrics_and_plots.R)
#   MAE  = sum|yhat - y| / N        (one division over all K gaps combined)
#   RMSE = sqrt( sum(yhat - y)^2 / N )
#   R²   = weighted mean of per-gap R²_k, weighted by n_k
#          (R² cannot be decomposed into additive sums across gaps)
#   where N = n_1 + n_2 + ... + n_K is the total number of half-hours across
#   all K artificial gaps in a given gap-size category.
#
# DEPENDENCIES
#   Requires compute_metrics_and_plots.R to have been run first so that
#   overall_metrics_pooled_NEE.csv exists under graphs/metrics_csv/ with
#   columns: site, management, model, model_key, target, gap_size,
#   mae_all, rmse_all, r2_all, mae_le30, rmse_le30, r2_le30,
#   mae_gt30, rmse_gt30, r2_gt30.
#
# SECTIONS
#   1  Settings        — paths, sites, models, variables, labels
#   2  Libraries       — packages and metric helper functions
#   3  Compute CSVs    — gap-level metrics -> Method I pooling -> write CSVs
#   4  Load BASE       — read BASE pooled metrics for delta computation
#   5  Compute deltas  — ΔMAE, ΔRMSE, ΔR², % MAE reduction (no-split only)
#   6  Heatmaps        — % MAE reduction tiles, one PNG per gap size
#   7  Summary         — output paths printed to console
# =============================================================================


# =============================================================================
# SECTION 1 — Settings
# =============================================================================

RESULTS_ROOT <- here::here("results")
GRAPHS_ROOT  <- here::here("graphs")
OUT_ROOT     <- file.path(GRAPHS_ROOT, "management_effect")

SITES       <- c("JC1", "JC2")
MODELS      <- c("RF", "MLP", "XGBoost")
TARGET_VARS <- c("NEE")
GAP_SIZES   <- c("S", "M", "L", "VL")
SPLIT_DAYS  <- 30   # kept for CSV completeness; not used in plots

MGMT_VARS <- c(
  "Grazing_days_since",
  "Fertiliser_days_since",
  "N",
  "grass_height",
  "grass_biomass",
  "phytomass_index"
)

# Symbol names from Table tab:variables — used as y-axis labels in heatmaps.
# Subscript notation (e.g. "dsg") is parsed by ggplot into proper math text.
MGMT_LABELS <- c(
  Grazing_days_since    = "dsg",
  Fertiliser_days_since = "dsf",
  N                     = "N",
  grass_height          = "h[s]",
  grass_biomass         = "B",
  phytomass_index       = "PI"
)

PRED_SUFFIX <- c(RF = "rf_predicted", MLP = "mlp_predicted", XGBoost = "xgb_predicted")


# =============================================================================
# SECTION 2 — Libraries & helpers
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(purrr)
  library(stringr); library(readr); library(ggplot2); library(tibble)
  library(forcats); library(lubridate)
})

make_dir <- function(...) dir.create(file.path(...), recursive = TRUE, showWarnings = FALSE)

# Method I metric helpers — identical to compute_metrics_and_plots.R
.mae  <- function(y, yh) { ok <- is.finite(y) & is.finite(yh); if (!any(ok)) return(NA_real_); mean(abs(yh[ok] - y[ok])) }
.rmse <- function(y, yh) { ok <- is.finite(y) & is.finite(yh); if (!any(ok)) return(NA_real_); sqrt(mean((yh[ok] - y[ok])^2)) }
.r2   <- function(y, yh) { ok <- is.finite(y) & is.finite(yh); if (sum(ok) < 2) return(NA_real_); cc <- suppressWarnings(cor(y[ok], yh[ok])); if (!is.finite(cc)) return(NA_real_); cc^2 }

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

MODEL_COLOURS <- c(RF = "#E41A1C", MLP = "#377EB8", XGBoost = "#4DAF4A")


# =============================================================================
# SECTION 3 — Compute & save management-effect CSVs
# =============================================================================
# For each site x mgmt_var x model x gap, compute raw error sums and per-gap
# R² for Method I pooled aggregation, then write one CSV per site x mgmt_var:
#   graphs/metrics_csv/management_effect/{site}/{mgmt_var}/
#     overall_metrics_pooled_NEE.csv

message("=== SECTION 3: Computing and saving management-effect CSVs ===")

compute_mgmt_gap_metrics <- function(site, model, mgmt_var) {
  rds_path <- file.path(
    RESULTS_ROOT, site, "management_effect", model, mgmt_var,
    "df_cv_all_predictions.rds"
  )
  if (!file.exists(rds_path)) { message("  Not found (skip): ", rds_path); return(tibble()) }
  
  preds <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.null(preds)) { message("  Load failed: ", rds_path); return(tibble()) }
  
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

message("  Computing gap-level metrics for all management runs ...")

all_mgmt_gap <- map_dfr(SITES, function(site) {
  map_dfr(MODELS, function(model) {
    map_dfr(MGMT_VARS, function(mv) {
      compute_mgmt_gap_metrics(site, model, mv)
    })
  })
}) %>%
  mutate(gap_size = factor(gap_size, levels = GAP_SIZES))

message("  Gap-level rows computed: ", nrow(all_mgmt_gap))

# Method I pooled aggregation
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

# Save CSVs
message("  Saving management-effect pooled CSVs ...")
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
    message("    Saved: ", site, " / ", mv)
  }
}
message("  All management-effect CSVs saved.")


# =============================================================================
# SECTION 4 — Load BASE metrics
# =============================================================================
# Reads graphs/metrics_csv/overall_metrics_pooled_NEE.csv
# Requires columns: site, management, model, target, gap_size, mae_all
# If site/management columns are missing, re-run compute_metrics_and_plots.R

message("\n=== SECTION 4: Loading BASE metrics ===")

base_csv_dir <- file.path(GRAPHS_ROOT, "metrics_csv")

base_pooled <- map_dfr(TARGET_VARS, function(tv) {
  p <- file.path(base_csv_dir, paste0("overall_metrics_pooled_", tv, ".csv"))
  if (!file.exists(p)) { message("  Not found: ", p); return(tibble()) }
  d <- read_csv(p, show_col_types = FALSE)
  # Guard: confirm required columns are present
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

message("  BASE rows loaded: ", nrow(base_pooled))


# =============================================================================
# SECTION 5 — Compute deltas (no-split only)
# =============================================================================
# delta_mae_all     = MAE(BASE+m) - MAE(BASE)          [negative = improvement]
# pct_mae_reduction = -100 * delta_mae_all / base_mae  [positive = improvement]

message("\n=== SECTION 5: Computing deltas ===")

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

message("  Delta rows: ", nrow(delta_pooled))


# =============================================================================
# SECTION 6 — Heatmaps (% MAE reduction, no split)
# =============================================================================
# Fill    : pct_mae_reduction = -100 * ΔMAE / MAE(BASE)
#           positive (green) = improvement; negative (red) = degradation
# Layout  : facet_wrap(~ site), y = management variable symbol, x = model
# y-axis  : symbol labels from Table tab:variables, parsed as math text
#           dsg,  dsf,  h[s] → h_s,  B, N, PI as-is
# Saved   : OUT_ROOT/heatmaps/NEE/pct_MAE_{gap_size}.png

message("\n=== SECTION 6: Heatmaps ===")

make_heatmap <- function(tv, gs) {
  d <- delta_pooled %>%
    filter(target == tv, gap_size == gs, !is.na(pct_mae_reduction)) %>%
    mutate(model = factor(model, levels = MODELS))
  if (!nrow(d)) return(NULL)
  
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
  message("  Heatmaps saved: ", tv)
}


# =============================================================================
# SECTION 7 — Summary
# =============================================================================

message("\n", strrep("=", 60))
message("COMPLETE")
message(strrep("=", 60))
message("CSVs  : graphs/metrics_csv/management_effect/{site}/{mgmt_var}/")
message("        overall_metrics_pooled_NEE.csv")
message("Plots : ", OUT_ROOT)
message("  1.  heatmaps/NEE/ — pct_MAE_{gap_size}.png")