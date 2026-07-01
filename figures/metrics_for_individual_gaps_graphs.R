# graph_colour_individual_gaps.R — Per-gap MAE/RMSE/R² plots, season & grazing encoded
#
# Reads the per-gap and pooled metric CSVs written by metrics_graphs.R
# (no recomputation) plus one reference df_cv_all_predictions.rds per site to recover
# each gap's dominant season and post-grazing recovery fraction. Draws two plot
# families across the four site × management panels: (A) split plots with ≤30 d / >30 d
# facet rows, and (B) whole-gap plots with no split. Line colour = dominant season,
# point shape = grazing-recovery fraction, black ◆ line = pooled (Method I) metric.
#
# Required input files (under graphs/metrics_csv/, from compute_metrics_and_plots.R):
#   gap_metrics_NEE_{S,M,L,VL}.csv        — per-gap points, ≤30 d / >30 d split
#   gap_metrics_whole_NEE_{S,M,L,VL}.csv  — per-gap points, whole gap (no split)
#   overall_metrics_pooled_NEE.csv        — pooled aggregate (black ◆ line)
#   one df_cv_all_predictions.rds per site — season (Spring/Summer/Autumn/Winter)
#                                            and Grazing_days_since metadata
#
# Outputs (PNG + multi-page PDF):
#   graphs/site_comparison_coloured/NEE/{S,M,L,VL}/  — split plots
#   graphs/site_comparison_whole/NEE/{S,M,L,VL}/     — whole-gap plots
#
# To adapt: change SITES_ALL / GAP_SIZES / TARGET_VARS, or the recovery thresholds
# (SPLIT_DAYS days, GRAZ_HIGH / GRAZ_MID fractions).


# Repository roots for reading metrics and writing graphs
RESULTS_ROOT <- here::here("results")
GRAPHS_ROOT  <- here::here("graphs")

# Post-grazing split (days) and the fraction cut-offs for the recovery shape groups
SPLIT_DAYS <- 30
GRAZ_HIGH  <- 0.80
GRAZ_MID   <- 0.50

# Target flux, gap-size classes, sites, and the fixed panel-column order
TARGET_VARS <- c("NEE")
GAP_SIZES   <- c("S", "M", "L", "VL")
SITES_ALL   <- c("JC1", "JC2")

COL_LABEL_LEVELS <- c("JC1 Unmanaged","JC1 Managed","JC2 Unmanaged","JC2 Managed")

# Resolve the df_cv_all_predictions.rds path for a given site / model / management
rds_path_for_model <- function(site, model_key, management) {
  if (model_key %in% c("minirec","mds")) {
    folder <- switch(model_key, minirec = "miniRECgap", mds = "MDS")
    file.path(RESULTS_ROOT, site, folder, "df_cv_all_predictions.rds")
  } else {
    folder <- switch(model_key, rf = "RF", mlp = "MLP", xgb = "XGBoost")
    file.path(RESULTS_ROOT, site, management, folder, "df_cv_all_predictions.rds")
  }
}


# Libraries
suppressPackageStartupMessages({
  library(here); library(dplyr); library(tibble); library(tidyr)
  library(purrr); library(stringr); library(readr); library(ggplot2)
})

# Create a (possibly nested) output directory if it does not exist
make_dir <- function(...) dir.create(file.path(...), recursive = TRUE, showWarnings = FALSE)

# Axis labels per metric and facet-row labels for the ≤30 d / >30 d windows
metric_ylab <- c(
  mae  = "MAE (\u03bcmol m\u207b\u00b2 s\u207b\u00b9)",
  rmse = "RMSE (\u03bcmol m\u207b\u00b2 s\u207b\u00b9)",
  r2   = "R\u00b2"
)
window_labels <- c(le30 = "\u226430 d since grazing", gt30 = ">30 d since grazing")

# Season colours and recovery-fraction point shapes
SEASON_COLOURS <- c(
  Spring       = "#4DAF4A",
  Summer       = "#FF7F00",
  Autumn       = "#A65628",
  Winter       = "#377EB8",
  Transitional = "#AAAAAA"
)
RECOVERY_SHAPES <- c(
  "High (>=80%)" = 21,
  "Mid (50-79%)" = 24,
  "Low (<50%)"   = 22
)

# Shared plot theme
base_theme <- theme_bw(base_size = 15) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.text        = element_text(size = 14),
    legend.position   = "bottom",
    legend.title      = element_text(size = 14),
    legend.text       = element_text(size = 14),
    plot.background   = element_rect(fill = "white", colour = NA),
    panel.spacing     = unit(1, "lines")
  )

# Which models appear in a panel: managed drops the benchmarks; MDS only at S/M
mods_for <- function(mgmt, gap_cat) {
  if (mgmt == "managed")         return(c("RF","MLP","XGBoost"))
  if (gap_cat %in% c("VL","L")) return(c("RF","MLP","XGBoost","miniRECgap"))
  c("RF","MLP","XGBoost","miniRECgap","MDS")
}


# Locate the metrics CSV folder produced by compute_metrics_and_plots.R
csv_dir <- file.path(GRAPHS_ROOT, "metrics_csv")
if (!dir.exists(csv_dir))
  stop("metrics_csv folder not found. Run compute_metrics_and_plots.R first.\nExpected: ", csv_dir)

# Per-gap dot values with the ≤30 d / >30 d split; MDS is dropped for L/VL
all_gap_metrics <- map_dfr(TARGET_VARS, function(tv) {
  map_dfr(GAP_SIZES, function(gs) {
    p <- file.path(csv_dir, paste0("gap_metrics_", tv, "_", gs, ".csv"))
    if (!file.exists(p)) { return(tibble()) }
    read_csv(p, show_col_types = FALSE)
  })
}) %>%
  mutate(
    model      = factor(model, levels = c("RF","MLP","XGBoost","miniRECgap","MDS")),
    gap_size   = factor(gap_size, levels = GAP_SIZES),
    management = factor(management, levels = c("unmanaged","managed"))
  ) %>%
  filter(!(model == "MDS" & gap_size %in% c("L","VL")))

# Per-gap dot values for the whole gap (no window split)
whole_gap_metrics <- map_dfr(TARGET_VARS, function(tv) {
  map_dfr(GAP_SIZES, function(gs) {
    p <- file.path(csv_dir, paste0("gap_metrics_whole_", tv, "_", gs, ".csv"))
    if (!file.exists(p)) { return(tibble()) }
    read_csv(p, show_col_types = FALSE)
  })
}) %>%
  mutate(
    model      = factor(model, levels = c("RF","MLP","XGBoost","miniRECgap","MDS")),
    gap_size   = factor(gap_size, levels = GAP_SIZES),
    management = factor(management, levels = c("unmanaged","managed"))
  ) %>%
  filter(!(model == "MDS" & gap_size %in% c("L","VL")))

# Pooled aggregate that becomes the black diamond line; columns prefixed agg_
agg_metrics <- map_dfr(TARGET_VARS, function(tv) {
  p <- file.path(csv_dir, paste0("overall_metrics_pooled_", tv, ".csv"))
  if (!file.exists(p)) { return(tibble()) }
  read_csv(p, show_col_types = FALSE)
}) %>%
  rename(
    agg_mae_all   = mae_all,
    agg_rmse_all  = rmse_all,
    agg_r2_all    = r2_all,
    agg_mae_le30  = mae_le30,
    agg_rmse_le30 = rmse_le30,
    agg_r2_le30   = r2_le30,
    agg_mae_gt30  = mae_gt30,
    agg_rmse_gt30 = rmse_gt30,
    agg_r2_gt30   = r2_gt30
  ) %>%
  mutate(
    model      = factor(model, levels = c("RF","MLP","XGBoost","miniRECgap","MDS")),
    gap_size   = factor(gap_size, levels = GAP_SIZES),
    management = factor(management, levels = c("unmanaged","managed"))
  ) %>%
  filter(!(model == "MDS" & gap_size %in% c("L","VL")))


# Pull per-gap season / grazing-recovery metadata from one reference RDS per site
extract_gap_meta_for_site <- function(site) {

  # Prefer RF managed, then RF/MLP unmanaged, then miniRECgap; take the first present
  candidates <- c(
    rds_path_for_model(site, "rf",      "managed"),
    rds_path_for_model(site, "rf",      "unmanaged"),
    rds_path_for_model(site, "mlp",     "unmanaged"),
    rds_path_for_model(site, "minirec", NA_character_)
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) {
    return(tibble())
  }

  preds <- tryCatch(readRDS(path), error = function(e) {
    NULL
  })
  if (is.null(preds)) return(tibble())

  # Gap-flag columns are named like S1, M3, L2, VL4
  flag_cols_all <- names(preds)[grepl("^(S|M|L|VL)\\d+$", names(preds))]
  if (!length(flag_cols_all)) {
    return(tibble())
  }

  has_season  <- all(c("Spring","Summer","Autumn","Winter") %in% names(preds))
  has_grazing <- "Grazing_days_since" %in% names(preds)

  # One metadata row per gap column
  map_dfr(flag_cols_all, function(fc) {
    gap_rows <- which(preds[[fc]] %in% c(TRUE, 1L))
    if (!length(gap_rows)) return(NULL)

    sub <- preds[gap_rows, ]

    # Dominant season = the one covering >50% of the gap's half-hours, else Transitional
    if (has_season) {
      props <- c(
        Spring = mean(sub$Spring == 1, na.rm = TRUE),
        Summer = mean(sub$Summer == 1, na.rm = TRUE),
        Autumn = mean(sub$Autumn == 1, na.rm = TRUE),
        Winter = mean(sub$Winter == 1, na.rm = TRUE)
      )
      dom <- if (max(props, na.rm = TRUE) > 0.5) names(which.max(props))
      else "Transitional"
    } else { dom <- "Transitional" }

    # Recovery group from the fraction of gap half-hours within SPLIT_DAYS of grazing
    if (has_grazing) {
      gd  <- as.numeric(sub$Grazing_days_since)
      pct <- mean(is.finite(gd) & gd <= SPLIT_DAYS, na.rm = TRUE)
    } else { pct <- 0 }
    rec <- dplyr::case_when(
      pct >= GRAZ_HIGH ~ "High (>=80%)",
      pct >= GRAZ_MID  ~ "Mid (50-79%)",
      TRUE             ~ "Low (<50%)"
    )

    gs_cat  <- str_extract(fc, "^[A-Z]+")
    gap_num <- as.integer(str_extract(fc, "\\d+$"))

    tibble(site = site, gap_size = gs_cat,
           gap_id = paste0(gs_cat, gap_num),
           dominant_season = dom, recovery_group = rec)
  })
}

# Metadata for both sites, with factor levels fixed for consistent colour/shape order
gap_meta <- map_dfr(SITES_ALL, extract_gap_meta_for_site) %>%
  mutate(
    gap_size        = factor(gap_size, levels = GAP_SIZES),
    dominant_season = factor(dominant_season,
                             levels = c("Spring","Summer","Autumn","Winter","Transitional")),
    recovery_group  = factor(recovery_group,
                             levels = c("High (>=80%)","Mid (50-79%)","Low (<50%)"))
  )

if (nrow(gap_meta) == 0)
  stop("No gap metadata extracted. Check that df_cv_all_predictions.rds files ",
       "contain gap flag columns and season/grazing columns.")


# Attach the season / recovery metadata to a metric table, filling missing as defaults
attach_meta <- function(df) {
  df %>%
    left_join(gap_meta %>% select(site, gap_size, gap_id,
                                  dominant_season, recovery_group),
              by = c("site", "gap_size", "gap_id")) %>%
    mutate(
      dominant_season = factor(
        replace_na(as.character(dominant_season), "Transitional"),
        levels = c("Spring","Summer","Autumn","Winter","Transitional")),
      recovery_group = factor(
        replace_na(as.character(recovery_group), "Low (<50%)"),
        levels = c("High (>=80%)","Mid (50-79%)","Low (<50%)"))
    )
}

all_gap_metrics   <- attach_meta(all_gap_metrics)
whole_gap_metrics <- attach_meta(whole_gap_metrics)


# Order models along x by the JC1-unmanaged ≤30 d reference metric (asc, or desc for R²)
model_order_global <- function(tv, gap_cat, metric) {
  ref_col <- paste0("agg_", metric, "_le30")
  agg_metrics %>%
    filter(target == tv, gap_size == gap_cat,
           site == "JC1", management == "unmanaged",
           !is.na(.data[[ref_col]])) %>%
    arrange(if (metric == "r2") desc(.data[[ref_col]]) else .data[[ref_col]]) %>%
    pull(model) %>% as.character()
}


# Build a split plot with ≤30 d / >30 d facet rows for one target/gap/metric
make_split_plot <- function(tv, gap_cat, metric,
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

  # Per-gap points: long over the two windows, one block per panel column
  d_indiv <- map_dfr(col_specs, function(spec) {
    mods <- mods_for(spec$mgmt, gap_cat)
    all_gap_metrics %>%
      filter(target == tv, gap_size == gap_cat,
             site == spec$site, management == spec$mgmt, model %in% mods) %>%
      select(model, gap_id, dominant_season, recovery_group,
             le30 = all_of(value_le), gt30 = all_of(value_gt)) %>%
      pivot_longer(c(le30, gt30), names_to = "window", values_to = "value") %>%
      mutate(col_label = paste(spec$site, str_to_title(spec$mgmt)))
  }) %>%
    filter(!is.na(value)) %>%
    mutate(
      window    = factor(window, levels = c("le30","gt30"), labels = window_labels),
      col_label = factor(col_label, levels = COL_LABEL_LEVELS),
      model     = factor(model, levels = ord_full)
    )

  # Matching pooled aggregate for the black diamond line
  d_agg <- map_dfr(col_specs, function(spec) {
    mods <- mods_for(spec$mgmt, gap_cat)
    agg_metrics %>%
      filter(target == tv, gap_size == gap_cat,
             site == spec$site, management == spec$mgmt, model %in% mods) %>%
      select(model, le30 = all_of(agg_le), gt30 = all_of(agg_gt)) %>%
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
  .build_plot(d_indiv, d_agg, metric, "window", COL_LABEL_LEVELS, add_title, title_str)
}


# Build a whole-gap plot (no window split) for one target/gap/metric
make_whole_plot <- function(tv, gap_cat, metric,
                            add_title = FALSE, title_str = "") {
  agg_col <- paste0("agg_", metric, "_all")

  ord_full <- c("RF", "XGBoost", "MLP", "miniRECgap", "MDS")

  col_specs <- list(
    list(site = "JC1", mgmt = "unmanaged"),
    list(site = "JC1", mgmt = "managed"),
    list(site = "JC2", mgmt = "unmanaged"),
    list(site = "JC2", mgmt = "managed")
  )

  # Per-gap points, one block per panel column
  d_indiv <- map_dfr(col_specs, function(spec) {
    mods <- mods_for(spec$mgmt, gap_cat)
    whole_gap_metrics %>%
      filter(target == tv, gap_size == gap_cat,
             site == spec$site, management == spec$mgmt, model %in% mods) %>%
      select(model, gap_id, dominant_season, recovery_group,
             value = all_of(metric)) %>%
      mutate(col_label = paste(spec$site, str_to_title(spec$mgmt)))
  }) %>%
    filter(!is.na(value)) %>%
    mutate(
      col_label = factor(col_label, levels = COL_LABEL_LEVELS),
      model     = factor(model, levels = ord_full)
    )

  # Matching pooled aggregate for the black diamond line
  d_agg <- map_dfr(col_specs, function(spec) {
    mods <- mods_for(spec$mgmt, gap_cat)
    agg_metrics %>%
      filter(target == tv, gap_size == gap_cat,
             site == spec$site, management == spec$mgmt, model %in% mods) %>%
      select(model, agg_value = all_of(agg_col)) %>%
      mutate(col_label = paste(spec$site, str_to_title(spec$mgmt)))
  }) %>%
    filter(!is.na(agg_value)) %>%
    mutate(
      col_label = factor(col_label, levels = COL_LABEL_LEVELS),
      model     = factor(model, levels = ord_full)
    )

  if (!nrow(d_indiv)) return(NULL)
  .build_plot(d_indiv, d_agg, metric, NULL, COL_LABEL_LEVELS, add_title, title_str)
}


# Assemble the ggplot shared by both families; facet_row_var = NULL means no split
.build_plot <- function(d_indiv, d_agg, metric, facet_row_var,
                        col_levels, add_title, title_str) {

  p <- ggplot(d_indiv, aes(x = model, y = value)) +
    geom_line(aes(group = gap_id, colour = dominant_season),
              linewidth = 0.45, alpha = 0.65) +
    geom_point(aes(fill = dominant_season, shape = recovery_group),
               colour = "black", stroke = 0.6,
               size  = 6.4,
               alpha = 0.3) +
    geom_line(data = d_agg,
              aes(x = model, y = agg_value, group = 1),
              colour = "black", linewidth = 1.1, inherit.aes = FALSE) +
    geom_point(data = d_agg,
               aes(x = model, y = agg_value),
               colour = "black", size = 6.5, shape = 18,
               inherit.aes = FALSE) +
    scale_colour_manual(name   = "Dominant season (>50%)",
                        values = SEASON_COLOURS, drop = FALSE,
                        breaks = c("Spring","Summer","Autumn","Winter")) +
    scale_fill_manual(name   = "Dominant season (>50%)",
                      values = SEASON_COLOURS, drop = FALSE,
                      breaks = c("Spring","Summer","Autumn","Winter")) +
    scale_shape_manual(name   = "Within 30 days of grazing",
                       values = RECOVERY_SHAPES, drop = FALSE) +
    scale_x_discrete(drop = TRUE) +
    labs(x = NULL, y = metric_ylab[[metric]],
         title = if (add_title) title_str else NULL) +
    base_theme +
    theme(
      axis.text.x  = element_text(angle = 35, hjust = 1, size = 14),
      strip.text.x = element_text(size = 14, face = "bold"),
      strip.text.y = element_text(angle = -90, size = 14, face = "bold")
    ) +
    guides(
      colour = guide_legend(
        title.position = "top",
        nrow = 1,
        override.aes = list(linewidth = 2, shape = NA)
      ),
      shape = guide_legend(title.position = "top", nrow = 1)
    )

  # space = "free_x" gives each column panel its own x-axis —
  # managed panels only show RF/MLP/XGBoost, miniRECgap is absent
  if (!is.null(facet_row_var)) {
    p <- p + facet_grid(reformulate("col_label", facet_row_var),
                        scales = "free", space = "free_x")
  } else {
    p <- p + facet_grid(. ~ col_label,
                        scales = "free", space = "free_x")
  }

  # R² gets 0.1 tick steps clamped to the data range; MAE/RMSE get integer 0.5 steps
  y_vals <- c(d_indiv$value, d_agg$agg_value)
  y_vals <- y_vals[is.finite(y_vals)]
  if (length(y_vals)) {
    if (metric == "r2") {
      y_lo <- floor(min(y_vals) * 10) / 10
      y_hi <- ceiling(max(y_vals) * 10) / 10
      p <- p + coord_cartesian(ylim = c(y_lo, y_hi)) +
        scale_y_continuous(breaks = seq(y_lo, y_hi, by = 0.1))
    } else {
      y_lo <- floor(min(y_vals))
      y_hi <- ceiling(max(y_vals))
      p <- p + scale_y_continuous(
        breaks       = seq(y_lo, y_hi, by = 0.5),
        labels       = function(x) ifelse(x == round(x), as.character(as.integer(x)), ""),
        minor_breaks = NULL
      )
    }
  }
  p
}


# Loop over targets and gap sizes, writing PNGs and one multi-page PDF per family
for (tv in TARGET_VARS) {
  for (gap_cat in GAP_SIZES) {

    # (A) Split plots (≤30 d / >30 d): one PNG per metric
    out_dir_split <- file.path(GRAPHS_ROOT, "site_comparison_coloured", tv, gap_cat)
    make_dir(out_dir_split)

    for (metric in c("mae","rmse","r2")) {
      p <- make_split_plot(tv, gap_cat, metric, add_title = FALSE)
      if (is.null(p)) next
      ggsave(file.path(out_dir_split, paste0(toupper(metric), "_", gap_cat, ".png")),
             p, width = 16, height = 10, dpi = 220, bg = "white")
    }

    # Split plots collected into a titled multi-page PDF
    pdf(file.path(out_dir_split, paste0("comparison_", gap_cat, "_", tv, ".pdf")),
        width = 16, height = 10, onefile = TRUE)
    for (metric in c("mae","rmse","r2")) {
      p <- make_split_plot(tv, gap_cat, metric, add_title = TRUE,
                           title_str = paste0(tv, " \u2014 ", gap_cat, " \u2014 ", toupper(metric),
                                              " | \u226430d / >30d split | colour = season",
                                              " | shape = grazing recovery | \u25c6 = pooled"))
      if (!is.null(p)) print(p)
    }
    dev.off()

    # (B) Whole-gap plots (no split): one PNG per metric
    out_dir_whole <- file.path(GRAPHS_ROOT, "site_comparison_whole", tv, gap_cat)
    make_dir(out_dir_whole)

    for (metric in c("mae","rmse","r2")) {
      p <- make_whole_plot(tv, gap_cat, metric, add_title = FALSE)
      if (is.null(p)) next
      ggsave(file.path(out_dir_whole, paste0(toupper(metric), "_", gap_cat, ".png")),
             p, width = 16, height = 6, dpi = 220, bg = "white")
    }

    # Whole-gap plots collected into a titled multi-page PDF
    pdf(file.path(out_dir_whole, paste0("comparison_whole_", gap_cat, "_", tv, ".pdf")),
        width = 16, height = 6, onefile = TRUE)
    for (metric in c("mae","rmse","r2")) {
      p <- make_whole_plot(tv, gap_cat, metric, add_title = TRUE,
                           title_str = paste0(tv, " \u2014 ", gap_cat, " \u2014 ", toupper(metric),
                                              " | whole gap (no split) | colour = season",
                                              " | shape = grazing recovery | \u25c6 = pooled"))
      if (!is.null(p)) print(p)
    }
    dev.off()
  }
}
