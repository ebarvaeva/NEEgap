# =============================================================================
# SITE COMPARISON PLOTS — Season & Recovery Encoded
# =============================================================================
#
# Reads CSVs produced by compute_metrics_and_plots.R — no recomputation:
#   gap_metrics_NEE_{gap_size}.csv        — per-gap dots (≤30d / >30d split)
#   gap_metrics_whole_NEE_{gap_size}.csv  — per-gap dots (whole gap, no split)
#   overall_metrics_pooled_NEE.csv        — black diamond aggregate line
#                                           (Method I pooled, one row per
#                                            site × management × model × gap_size)
#
# Then reads one reference df_cv_all_predictions.rds per site to extract
# per-gap season and grazing-recovery metadata (shape / colour encoding).
#
# LINE COLOUR — dominant season of the gap
#   Spring (green), Summer (orange), Autumn (brown), Winter (blue),
#   Transitional (grey) when no season > 50 % of gap half-hours.
#
# DOT SHAPE — fraction of gap half-hours within 30 d of last grazing
#   Circle (●) >= 80 %  |  Triangle (▲) 50–79 %  |  Square (■) < 50 %
#
# Black ◆ aggregate line = Method I pooled metric from overall_metrics_pooled_NEE.csv
#
# TWO PLOT FAMILIES:
#   (A) Split plots   — facet rows: ≤30d / >30d since grazing
#                       OUTPUT: graphs/site_comparison_coloured/NEE/{GAP_SIZE}/
#   (B) Whole-gap plots — no row split; one panel row per column
#                         OUTPUT: graphs/site_comparison_whole/NEE/{GAP_SIZE}/
#
# =============================================================================


# =============================================================================
# SECTION 1 — Settings
# =============================================================================

RESULTS_ROOT <- here::here("results")
GRAPHS_ROOT  <- here::here("graphs")

SPLIT_DAYS <- 30
GRAZ_HIGH  <- 0.80
GRAZ_MID   <- 0.50

TARGET_VARS <- c("NEE")
GAP_SIZES   <- c("S", "M", "L", "VL")
SITES_ALL   <- c("JC1", "JC2")

COL_LABEL_LEVELS <- c("JC1 Unmanaged","JC1 Managed","JC2 Unmanaged","JC2 Managed")

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

mods_for <- function(mgmt, gap_cat) {
  if (mgmt == "managed")         return(c("RF","MLP","XGBoost"))
  if (gap_cat %in% c("VL","L")) return(c("RF","MLP","XGBoost","miniRECgap"))
  c("RF","MLP","XGBoost","miniRECgap","MDS")
}


# =============================================================================
# SECTION 3 — Load CSVs from compute_metrics_and_plots.R
# =============================================================================

message("Reading gap metrics CSVs ...")
csv_dir <- file.path(GRAPHS_ROOT, "metrics_csv")
if (!dir.exists(csv_dir))
  stop("metrics_csv folder not found. Run compute_metrics_and_plots.R first.\nExpected: ", csv_dir)

# ---- Per-gap dot values — split (≤30d / >30d) ----
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
  filter(!(model == "MDS" & gap_size %in% c("L","VL")))

message("  Loaded ", nrow(all_gap_metrics), " per-gap split metric rows.")

# ---- Per-gap dot values — whole gap (no split) ----
whole_gap_metrics <- map_dfr(TARGET_VARS, function(tv) {
  map_dfr(GAP_SIZES, function(gs) {
    p <- file.path(csv_dir, paste0("gap_metrics_whole_", tv, "_", gs, ".csv"))
    if (!file.exists(p)) { message("  Missing: ", p); return(tibble()) }
    read_csv(p, show_col_types = FALSE)
  })
}) %>%
  mutate(
    model      = factor(model, levels = c("RF","MLP","XGBoost","miniRECgap","MDS")),
    gap_size   = factor(gap_size, levels = GAP_SIZES),
    management = factor(management, levels = c("unmanaged","managed"))
  ) %>%
  filter(!(model == "MDS" & gap_size %in% c("L","VL")))

message("  Loaded ", nrow(whole_gap_metrics), " per-gap whole-gap metric rows.")

# ---- Pooled aggregate for black diamond line ----
agg_metrics <- map_dfr(TARGET_VARS, function(tv) {
  p <- file.path(csv_dir, paste0("overall_metrics_pooled_", tv, ".csv"))
  if (!file.exists(p)) { message("  Missing: ", p); return(tibble()) }
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

message("  Loaded ", nrow(agg_metrics), " pooled aggregate rows (black diamond line).")


# =============================================================================
# SECTION 4 — Extract gap metadata from one reference RDS per site
# =============================================================================

message("Extracting gap metadata (season / grazing recovery) ...")

extract_gap_meta_for_site <- function(site) {
  
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
    message("  Season columns absent — dominant_season = Transitional")
  if (!has_grazing)
    message("  Grazing_days_since absent — recovery_group = Low (<50%)")
  
  map_dfr(flag_cols_all, function(fc) {
    gap_rows <- which(preds[[fc]] %in% c(TRUE, 1L))
    if (!length(gap_rows)) return(NULL)
    
    sub <- preds[gap_rows, ]
    
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

gap_meta <- map_dfr(SITES_ALL, extract_gap_meta_for_site) %>%
  mutate(
    gap_size        = factor(gap_size, levels = GAP_SIZES),
    dominant_season = factor(dominant_season,
                             levels = c("Spring","Summer","Autumn","Winter","Transitional")),
    recovery_group  = factor(recovery_group,
                             levels = c("High (>=80%)","Mid (50-79%)","Low (<50%)"))
  )

message("  Metadata rows: ", nrow(gap_meta),
        " (", n_distinct(gap_meta$gap_id), " unique gap IDs per site)")

if (nrow(gap_meta) == 0)
  stop("No gap metadata extracted. Check that df_cv_all_predictions.rds files ",
       "contain gap flag columns and season/grazing columns.")


# =============================================================================
# SECTION 5 — Join metadata to both metric tables
# =============================================================================

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

message("  Join complete. Split rows: ", nrow(all_gap_metrics),
        "  Whole-gap rows: ", nrow(whole_gap_metrics))


# =============================================================================
# SECTION 6 — Model order helper
# =============================================================================

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
# SECTION 7A — Plot function: split (≤30d / >30d rows)
# =============================================================================

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


# =============================================================================
# SECTION 7B — Plot function: whole gap (no 30d split)
# =============================================================================

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


# =============================================================================
# SECTION 7C — Shared plot builder
# =============================================================================

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


# =============================================================================
# SECTION 8 — Generate plots
# =============================================================================

message("\nGenerating plots ...")

for (tv in TARGET_VARS) {
  for (gap_cat in GAP_SIZES) {
    
    # ---- (A) Split plots (≤30d / >30d) --------------------------------------
    out_dir_split <- file.path(GRAPHS_ROOT, "site_comparison_coloured", tv, gap_cat)
    make_dir(out_dir_split)
    
    for (metric in c("mae","rmse","r2")) {
      p <- make_split_plot(tv, gap_cat, metric, add_title = FALSE)
      if (is.null(p)) next
      ggsave(file.path(out_dir_split, paste0(toupper(metric), "_", gap_cat, ".png")),
             p, width = 16, height = 10, dpi = 220, bg = "white")
      message("  Saved (split): ", toupper(metric), "_", gap_cat, ".png")
    }
    
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
    message("  Saved PDF (split): comparison_", gap_cat, "_", tv, ".pdf")
    
    # ---- (B) Whole-gap plots (no split) -------------------------------------
    out_dir_whole <- file.path(GRAPHS_ROOT, "site_comparison_whole", tv, gap_cat)
    make_dir(out_dir_whole)
    
    for (metric in c("mae","rmse","r2")) {
      p <- make_whole_plot(tv, gap_cat, metric, add_title = FALSE)
      if (is.null(p)) next
      ggsave(file.path(out_dir_whole, paste0(toupper(metric), "_", gap_cat, ".png")),
             p, width = 16, height = 6, dpi = 220, bg = "white")
      message("  Saved (whole): ", toupper(metric), "_", gap_cat, ".png")
    }
    
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
    message("  Saved PDF (whole): comparison_whole_", gap_cat, "_", tv, ".pdf")
  }
}

message("\n=== Done.")
message("Split plots : ", file.path(GRAPHS_ROOT, "site_comparison_coloured"))
message("Whole plots : ", file.path(GRAPHS_ROOT, "site_comparison_whole"))