# =============================================================================
# graph_colour_individual_gaps.R — Season and Recovery Encoded Comparison Plots
# =============================================================================
#
# PURPOSE
#   Extends the connected-dot plots from compute_metrics_and_plots.R by
#   encoding two additional dimensions of gap metadata into the visual:
#
#   LINE COLOUR — dominant season of the gap
#     Each gap label is assigned the season that accounts for more than 50%
#     of its half-hourly timestamps.  If no single season exceeds 50%, the
#     gap is labelled "Transitional".
#       Spring = green   Summer = orange   Autumn = brown
#       Winter = blue    Transitional = grey
#
#   DOT SHAPE — grazing recovery phase of the gap
#     Quantified as the fraction of gap half-hours within SPLIT_DAYS of the
#     most recent grazing event (Grazing_days_since <= SPLIT_DAYS).
#       Circle  ● — High  (>= 80%  of half-hours in early recovery)
#       Triangle▲ — Mid   (50–79% in early recovery)
#       Square  ■ — Low   (< 50%  in early recovery)
#
#   These encodings allow the reader to identify whether certain seasons or
#   recovery phases drive systematic differences between models, which is
#   not visible in the plain grey-line version from compute_metrics_and_plots.R.
#
# DESIGN PRINCIPLE: METADATA IS A GAP PROPERTY, NOT A MODEL PROPERTY
#   Season and recovery phase are intrinsic properties of each artificial gap
#   (they depend only on the timestamps and Grazing_days_since in the raw data,
#   not on model predictions).  A single reference RDS per site is therefore
#   sufficient to extract the metadata; it is then joined on site × gap_id
#   and shared across all models.  This guarantees that the same gap label
#   carries the same visual encoding in all four column panels.
#
# INPUTS
#   graphs/metrics_csv/gap_metrics_{TARGET}_{GAP_SIZE}.csv
#     Must exist — run compute_metrics_and_plots.R first.
#   results/{SITE}/{MANAGEMENT}/{MODEL}/df_cv_all_predictions.rds
#     One reference RDS per site used to extract gap metadata.
#
# OUTPUTS
#   graphs/site_comparison_coloured/{TARGET}/{GAP_SIZE}/
#     {MAE|RMSE|R2}_{GAP_SIZE}.png   — PNG (no title)
#     comparison_{GAP_SIZE}_{TARGET}.pdf — PDF (with title)
#
# SECTIONS
#   1  Settings          — thresholds, target vars, gap sizes, column layout
#   2  Libraries         — packages, theme, colour/shape palettes
#   3  Load CSV metrics  — read gap_metrics CSVs; compute weighted aggregate
#   4  Extract gap metadata — read reference RDS per site; compute season/recovery
#   5  Join metadata     — attach season/recovery to all_gap_metrics by gap_id
#   6  Model order       — best-to-worst ordering from JC1 unmanaged aggregate
#   7  Plot function     — make_encoded_comparison_plot() definition
#   8  Generate plots    — save loop over target × gap size × metric
#
# HOW TO RUN
#   source("metrics/graph_colour_individual_gaps.R")
#   compute_metrics_and_plots.R must have been run first.
#
# =============================================================================

# =============================================================================
# SECTION 1 — Settings
# =============================================================================
# Edit these variables to adapt the script to a different analysis.

RESULTS_ROOT <- here::here("results")
GRAPHS_ROOT  <- here::here("graphs")

# Grazing recovery thresholds (fraction of gap half-hours within SPLIT_DAYS)
SPLIT_DAYS <- 30    # days since grazing defining early vs late recovery
GRAZ_HIGH  <- 0.80  # >= this fraction → "High" recovery group (circle dots)
GRAZ_MID   <- 0.50  # >= this fraction → "Mid"  recovery group (triangle dots)
                    # < GRAZ_MID        → "Low"  recovery group (square dots)

TARGET_VARS <- c("NEE", "Reco", "GPP")
GAP_SIZES   <- c("S", "M", "L", "VL")
SITES_ALL   <- c("JC1", "JC2")

# Column panel order — must match COL_LABEL_LEVELS in compute_metrics_and_plots.R
COL_LABEL_LEVELS <- c("JC1 Unmanaged","JC1 Managed","JC2 Unmanaged","JC2 Managed")

# Model catalogue — must be identical to MODEL_CATALOG in compute_metrics_and_plots.R
# so that model labels and management assignments are consistent across scripts.
MODEL_CATALOG <- tibble::tribble(
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

rds_path_for_model <- function(site, model_key, management) {
  if (model_key %in% c("minirec","mds")) {
    folder <- switch(model_key, minirec = "miniRECgap", mds = "MDS")
    file.path(RESULTS_ROOT, site, folder, "df_cv_all_predictions.rds")
  } else {
    folder <- switch(model_key, rf = "RF", mlp = "MLP", xgb = "XGBoost")
    file.path(RESULTS_ROOT, site, management, folder, "df_cv_all_predictions.rds")
  }
}

# =============================================================================
# SECTION 2 — Libraries
# =============================================================================
# Colour and shape palettes are defined here as named vectors so that the
# same visual encodings apply regardless of which seasons or recovery groups
# happen to be present in any particular subset of data.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tibble); library(tidyr)
  library(purrr); library(stringr); library(readr); library(ggplot2)
})

make_dir <- function(...) dir.create(file.path(...), recursive = TRUE, showWarnings = FALSE)

metric_ylab <- c(
  mae  = "MAE (\u03bcmol m\u207b\u00b2 s\u207b\u00b9)",
  rmse = "RMSE (\u03bcmol m\u207b\u00b2 s\u207b\u00b9)",
  r2   = "R\u00b2"
)
window_labels <- c(le30 = "\u226430 d since grazing", gt30 = ">30 d since grazing")

# Season colour palette — colour-blind friendly set matched to seasonal ecology
# (green for growth, orange for harvest, brown for senescence, blue for dormancy)
SEASON_COLOURS <- c(
  Spring       = "#4DAF4A",
  Summer       = "#FF7F00",
  Autumn       = "#A65628",
  Winter       = "#377EB8",
  Transitional = "#AAAAAA"
)
# Recovery phase shapes:  circle = high recovery intensity,
# triangle = medium, square = low (echoes the gradient from intense to low impact)
RECOVERY_SHAPES <- c(
  "High (>=80%)" = 16,
  "Mid (50-79%)" = 17,
  "Low (<50%)"   = 15
)

base_theme <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(size = 10),
    legend.position   = "bottom",
    legend.title      = element_text(size = 9),
    legend.text       = element_text(size = 8),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.spacing     = unit(1, "lines")
  )

# Returns the set of models valid for a given management condition and gap size.
# MDS is excluded from L and VL gaps (unreliable for long gaps).
# miniRECgap and MDS have no managed-condition run.
mods_for <- function(mgmt, gap_cat) {
  if (mgmt == "managed")           return(c("RF","MLP","XGBoost"))
  if (gap_cat %in% c("VL","L"))   return(c("RF","MLP","XGBoost","miniRECgap"))
  c("RF","MLP","XGBoost","miniRECgap","MDS")
}

# =============================================================================
# SECTION 3 — Load CSV metrics
# =============================================================================
# Reads the per-gap metric CSVs produced by compute_metrics_and_plots.R.
# No metric recomputation is performed here; the CSVs are the canonical source.
# The weighted aggregate (agg_metrics) is recomputed from these CSVs for the
# black reference line in the plots.

message("Reading gap metrics CSVs ...")
csv_dir <- file.path(GRAPHS_ROOT, "metrics_csv")
if (!dir.exists(csv_dir))
  stop("metrics_csv folder not found. Run compute_metrics_and_plots.R first.\nExpected: ", csv_dir)

all_gap_metrics <- map_dfr(TARGET_VARS, function(tv) {
  map_dfr(GAP_SIZES, function(gs) {
    p <- file.path(csv_dir, paste0("gap_metrics_", tv, "_", gs, ".csv"))
    if (!file.exists(p)) { message("  Missing: ", p); return(tibble()) }
    read_csv(p, show_col_types = FALSE)
  })
}) %>%
  mutate(
    model      = factor(model, levels = c("RF","MLP","XGBoost","miniRECgap","MDS")),
    gap_size   = factor(gap_size, levels = GAP_SIZES),
    management = factor(management, levels = c("unmanaged","managed"))
  ) %>%
  filter(!(model == "MDS" & gap_size == "VL"))

message("  Loaded ", nrow(all_gap_metrics), " metric rows.")

# Aggregate for black reference line
agg_metrics <- all_gap_metrics %>%
  group_by(site, management, model, target, gap_size) %>%
  summarise(
    agg_mae_le30  = weighted.mean(mae_le30,  n_le30, na.rm = TRUE),
    agg_rmse_le30 = weighted.mean(rmse_le30, n_le30, na.rm = TRUE),
    agg_r2_le30   = weighted.mean(r2_le30,   n_le30, na.rm = TRUE),
    agg_mae_gt30  = weighted.mean(mae_gt30,  n_gt30, na.rm = TRUE),
    agg_rmse_gt30 = weighted.mean(rmse_gt30, n_gt30, na.rm = TRUE),
    agg_r2_gt30   = weighted.mean(r2_gt30,   n_gt30, na.rm = TRUE),
    .groups = "drop"
  )

# =============================================================================
# SECTION 4 — Extract gap metadata from one reference RDS per site
# =============================================================================
# Season and recovery phase are properties of the gap (determined by its
# timestamps and Grazing_days_since), not of any particular model.  A single
# reference RDS per site is therefore sufficient.
#
# extract_gap_meta_for_site() tries managed RF → unmanaged RF → unmanaged MLP
# → miniRECgap in order and uses the first file that exists.  For each gap
# flag column, it computes:
#   dominant_season  — season exceeding 50% of gap half-hours, or "Transitional"
#   recovery_group   — "High" / "Mid" / "Low" based on Grazing_days_since fraction
#
# The metadata is then joined onto all_gap_metrics by site + gap_size + gap_id
# in Section 5 so that all models share the same encoding for the same gap.

message("Extracting gap metadata (one reference RDS per site) ...")

extract_gap_meta_for_site <- function(site) {
  
  # Try managed RF first, then unmanaged RF, then any available RDS
  candidates <- c(
    rds_path_for_model(site, "rf",      "managed"),
    rds_path_for_model(site, "rf",      "unmanaged"),
    rds_path_for_model(site, "mlp",     "unmanaged"),
    rds_path_for_model(site, "minirec", NA_character_)
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) {
    message("  No RDS found for site ", site, " — gap metadata unavailable.")
    return(tibble())
  }
  message("  Using reference RDS for ", site, ": ", basename(dirname(path)))
  
  preds <- tryCatch(readRDS(path), error = function(e) {
    message("  Failed: ", path, " — ", conditionMessage(e)); NULL
  })
  if (is.null(preds)) return(tibble())
  
  flag_cols_all <- names(preds)[grepl("^(S|M|L|VL)\\d+$", names(preds))]
  if (!length(flag_cols_all)) {
    message("  No gap flag columns in: ", path); return(tibble())
  }
  
  has_season  <- all(c("Spring","Summer","Autumn","Winter") %in% names(preds))
  has_grazing <- "Grazing_days_since" %in% names(preds)
  
  if (!has_season)
    message("  Season columns absent in ", path, " — dominant_season = Transitional")
  if (!has_grazing)
    message("  Grazing_days_since absent in ", path, " — recovery_group = Low (<50%)")
  
  map_dfr(flag_cols_all, function(fc) {
    gap_rows <- which(preds[[fc]] %in% c(TRUE, 1L))
    if (!length(gap_rows)) return(NULL)
    
    sub <- preds[gap_rows, ]
    
    # Dominant season — computed over ALL half-hours of the gap
    if (has_season) {
      props <- c(
        Spring = mean(sub$Spring == 1, na.rm = TRUE),
        Summer = mean(sub$Summer == 1, na.rm = TRUE),
        Autumn = mean(sub$Autumn == 1, na.rm = TRUE),
        Winter = mean(sub$Winter == 1, na.rm = TRUE)
      )
      dom <- if (max(props, na.rm = TRUE) > 0.5) names(which.max(props))
      else "Transitional"
    } else {
      dom <- "Transitional"
    }
    
    # Recovery group — % of ALL gap half-hours within SPLIT_DAYS of last grazing
    # (not split by le30/gt30 — the shape is the same for both panel rows)
    if (has_grazing) {
      gd  <- as.numeric(sub$Grazing_days_since)
      pct <- mean(is.finite(gd) & gd <= SPLIT_DAYS, na.rm = TRUE)
    } else {
      pct <- 0
    }
    rec <- dplyr::case_when(
      pct >= GRAZ_HIGH ~ "High (>=80%)",
      pct >= GRAZ_MID  ~ "Mid (50-79%)",
      TRUE             ~ "Low (<50%)"
    )
    
    gs_cat  <- str_extract(fc, "^[A-Z]+")
    gap_num <- as.integer(str_extract(fc, "\\d+$"))
    
    tibble(
      site            = site,
      gap_size        = gs_cat,
      gap_id          = paste0(gs_cat, gap_num),
      dominant_season = dom,
      recovery_group  = rec
    )
  })
}

# One call per site — metadata is site-level, not model-level
gap_meta <- map_dfr(SITES_ALL, extract_gap_meta_for_site) %>%
  mutate(
    gap_size        = factor(gap_size, levels = GAP_SIZES),
    dominant_season = factor(dominant_season,
                             levels = c("Spring","Summer","Autumn","Winter","Transitional")),
    recovery_group  = factor(recovery_group,
                             levels = c("High (>=80%)","Mid (50-79%)","Low (<50%)"))
  )

message("  Metadata rows: ", nrow(gap_meta), " (", n_distinct(gap_meta$gap_id), " unique gap IDs per site)")

if (nrow(gap_meta) == 0)
  stop("No gap metadata extracted. Check that df_cv_all_predictions.rds files ",
       "contain gap flag columns (S1, M2, ...) and season/grazing columns.")

# =============================================================================
# SECTION 5 — Join metadata onto metrics
# =============================================================================
# Left-join gap_meta onto all_gap_metrics by site × gap_size × gap_id.
# The join key is the gap identity, not the model, ensuring that all rows for
# the same gap carry the same season and recovery group regardless of model.
# Unmatched rows (join misses, or gaps in unmanaged sites without grazing data)
# default to "Transitional" season and "Low (<50%)" recovery group.

all_gap_metrics <- all_gap_metrics %>%
  left_join(
    gap_meta %>% select(site, gap_size, gap_id, dominant_season, recovery_group),
    by = c("site", "gap_size", "gap_id")
  ) %>%
  mutate(
    dominant_season = factor(
      replace_na(as.character(dominant_season), "Transitional"),
      levels = c("Spring","Summer","Autumn","Winter","Transitional")
    ),
    recovery_group = factor(
      replace_na(as.character(recovery_group), "Low (<50%)"),
      levels = c("High (>=80%)","Mid (50-79%)","Low (<50%)")
    )
  )

n_no_meta <- sum(all_gap_metrics$dominant_season == "Transitional" &
                   all_gap_metrics$recovery_group == "Low (<50%)")
message("  Rows defaulting to Transitional/Low (join miss or genuine): ",
        n_no_meta, " of ", nrow(all_gap_metrics))

# =============================================================================
# SECTION 6 — Model order helper
# =============================================================================
# model_order_global() ranks models from best to worst for a given target,
# gap size, and metric using the JC1 unmanaged aggregate as the reference.
# This ordering is applied to the x-axis of all panels so that the best model
# is consistently leftmost across all four column panels.

model_order_global <- function(tv, gap_cat, metric) {
  ref_col <- paste0("agg_", metric, "_le30")
  agg_metrics %>%
    filter(target == tv, gap_size == gap_cat,
           site == "JC1", management == "unmanaged",
           !is.na(.data[[ref_col]])) %>%
    arrange(if (metric == "r2") desc(.data[[ref_col]]) else .data[[ref_col]]) %>%
    pull(model) %>% as.character()
}

# =============================================================================
# SECTION 7 — Plot function — make_encoded_comparison_plot()
# =============================================================================
# Extends the plain connected-dot plot from compute_metrics_and_plots.R with:
#   - Line colour = dominant_season (one colour per season, grey if transitional)
#   - Dot shape = recovery_group (circle/triangle/square)
# The black diamond aggregate line and the grey individual-gap line structure
# are retained from the original.  The y-axis range is computed from the data
# and rounded to neat boundaries (0.1 for R², 1 unit for MAE/RMSE).
# The function returns NULL silently if no data is available, allowing the
# save loop to skip missing combinations without error.

make_encoded_comparison_plot <- function(tv, gap_cat, metric,
                                         add_title = FALSE, title_str = "") {
  value_le <- paste0(metric, "_le30")
  value_gt <- paste0(metric, "_gt30")
  agg_le   <- paste0("agg_", metric, "_le30")
  agg_gt   <- paste0("agg_", metric, "_gt30")
  
  ord      <- model_order_global(tv, gap_cat, metric)
  all_mods <- unique(c(mods_for("unmanaged", gap_cat), mods_for("managed", gap_cat)))
  ord_full <- c(ord, setdiff(all_mods, ord))
  
  col_specs <- list(
    list(site = "JC1", mgmt = "unmanaged"),
    list(site = "JC1", mgmt = "managed"),
    list(site = "JC2", mgmt = "unmanaged"),
    list(site = "JC2", mgmt = "managed")
  )
  
  d_indiv <- map_dfr(col_specs, function(spec) {
    mods <- mods_for(spec$mgmt, gap_cat)
    all_gap_metrics %>%
      filter(target == tv, gap_size == gap_cat,
             site == spec$site,
             (management == spec$mgmt | (is.na(management) & spec$mgmt == "unmanaged")),
             model %in% mods) %>%
      select(model, gap_id, dominant_season, recovery_group,
             le30 = all_of(value_le),
             gt30 = all_of(value_gt)) %>%
      pivot_longer(c(le30, gt30), names_to = "window", values_to = "value") %>%
      mutate(col_label = paste(spec$site, str_to_title(spec$mgmt)))
  }) %>%
    filter(!is.na(value)) %>%
    mutate(
      window    = factor(window, levels = c("le30","gt30"), labels = window_labels),
      col_label = factor(col_label, levels = COL_LABEL_LEVELS),
      model     = factor(model, levels = ord_full)
    )
  
  d_agg <- map_dfr(col_specs, function(spec) {
    mods <- mods_for(spec$mgmt, gap_cat)
    agg_metrics %>%
      filter(target == tv, gap_size == gap_cat,
             site == spec$site, management == spec$mgmt,
             model %in% mods) %>%
      select(model,
             le30 = all_of(agg_le),
             gt30 = all_of(agg_gt)) %>%
      pivot_longer(c(le30, gt30), names_to = "window", values_to = "agg_value") %>%
      mutate(col_label = paste(spec$site, str_to_title(spec$mgmt)))
  }) %>%
    filter(!is.na(agg_value)) %>%
    mutate(
      window    = factor(window, levels = c("le30","gt30"), labels = window_labels),
      col_label = factor(col_label, levels = COL_LABEL_LEVELS),
      model     = factor(model, levels = ord_full)
    )
  
  if (!nrow(d_indiv)) return(NULL)
  
  p <- ggplot(d_indiv, aes(x = model, y = value)) +
    geom_line(aes(group = gap_id, colour = dominant_season),
              linewidth = 0.45, alpha = 0.65) +
    geom_point(aes(colour = dominant_season, shape = recovery_group),
               size = 2.2, alpha = 0.9) +
    geom_line(data = d_agg,
              aes(x = model, y = agg_value, group = 1),
              colour = "black", linewidth = 1.1, inherit.aes = FALSE) +
    geom_point(data = d_agg,
               aes(x = model, y = agg_value),
               colour = "black", size = 3.5, shape = 18, inherit.aes = FALSE) +
    facet_grid(window ~ col_label, scales = "free") +
    scale_colour_manual(name   = "Dominant season (>50%)",
                        values = SEASON_COLOURS, drop = FALSE,
                        breaks = c("Spring","Summer","Autumn","Winter")) +
    scale_shape_manual(name    = "Within 30 days of grazing",
                       values  = RECOVERY_SHAPES, drop = FALSE) +
    scale_x_discrete(drop = TRUE) +
    labs(x = NULL, y = metric_ylab[[metric]],
         title = if (add_title) title_str else NULL) +
    base_theme +
    theme(
      axis.text.x  = element_text(angle = 35, hjust = 1, size = 8),
      strip.text.x = element_text(size = 9, face = "bold"),
      strip.text.y = element_text(angle = -90, size = 9, face = "bold")
    ) +
    guides(
      colour = guide_legend(title.position = "top", nrow = 1,
                            override.aes = list(linewidth = 2, shape = NA)),
      shape  = guide_legend(title.position = "top", nrow = 1)
    )
  
  y_vals <- c(d_indiv$value, d_agg$agg_value)
  y_vals <- y_vals[is.finite(y_vals)]
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

# =============================================================================
# SECTION 8 — Generate plots
# =============================================================================
# Loops over all target × gap size combinations and saves both PNG and PDF.
# PNG has no title (clean for paper figures); PDF has a full descriptive title
# including the visual encoding legend for review purposes.

message("\nGenerating encoded site-comparison plots ...")

for (tv in TARGET_VARS) {
  for (gap_cat in GAP_SIZES) {
    
    out_dir <- file.path(GRAPHS_ROOT, "site_comparison_coloured", tv, gap_cat)
    make_dir(out_dir)
    
    for (metric in c("mae","rmse","r2")) {
      p <- make_encoded_comparison_plot(tv, gap_cat, metric, add_title = FALSE)
      if (is.null(p)) next
      out_png <- file.path(out_dir, paste0(toupper(metric), "_", gap_cat, ".png"))
      ggsave(out_png, p, width = 16, height = 7, dpi = 220, bg = "white")
      message("  Saved: ", out_png)
    }
    
    out_pdf <- file.path(out_dir, paste0("comparison_", gap_cat, "_", tv, ".pdf"))
    pdf(out_pdf, width = 16, height = 7, onefile = TRUE)
    for (metric in c("mae","rmse","r2")) {
      p <- make_encoded_comparison_plot(
        tv, gap_cat, metric, add_title = TRUE,
        title_str = paste0(tv, " \u2014 ", gap_cat, " \u2014 ", toupper(metric),
                           " | colour=season | shape=grazing recovery | \u25c6=overall")
      )
      if (!is.null(p)) print(p)
    }
    dev.off()
    message("  Saved PDF: ", out_pdf)
  }
}

message("\n=== Done. Output: ", file.path(GRAPHS_ROOT, "site_comparison_coloured"))

