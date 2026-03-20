# =============================================================================
# variable_importance_graphs.R — RF Variable Importance Plots
# =============================================================================
#
# PURPOSE
#   Produces all plots needed to evaluate the relative importance of
#   management variables as RF predictors across managed/unmanaged sites
#   and the management variable ablation study.
#
# BACKGROUND
#   Random Forest (ranger) saves impurity-based variable importance
#   (mean decrease in node impurity, Gini criterion) for each predictor
#   at every split across all trees.  For each artificial gap label,
#   the raw importance vector is normalised to relative importance
#   (importance_rel = importance / sum(importance)), giving the fraction
#   of total impurity attributable to each predictor.  These values are
#   averaged across gap labels to produce the mean relative importance
#   displayed in all plots.
#
# THREE DATA FAMILIES
# ───────────────────
#  A. MANAGED
#     results/{SITE}/managed/RF/variable_importance/
#     Importance under the full managed predictor set (BASE + all management).
#
#  B. UNMANAGED
#     results/{SITE}/unmanaged/RF/variable_importance/
#     Importance under the BASE-only predictor set.  Serves as the reference
#     showing which meteorological and temporal predictors dominate when no
#     management information is available.
#
#  C. MANAGEMENT EFFECT
#     results/{SITE}/management_effect/RF/{MGMT_VAR}/variable_importance/
#     Importance when exactly one management variable is added to BASE.
#     Shows where each management variable ranks relative to meteorological
#     predictors and how it redistributes importance among them.
#     The Phytomass Index (PI) run is included under phytomass_index/.
#
# NOTE: Variable importance is saved only by RF_CV.R.  Figures 3, 4 (RF only)
# and 5 (RF, gap size M) use per-gap-label detail from the .rds files.
# Figures A1, A2, 1, 2 are based on the aggregate _all.rds files and are
# model-agnostic in principle, but currently only RF importance is loaded.
#
# FIGURES PRODUCED
# ────────────────
#  Family A/B  →  graphs/variable_importance/00_managed/  00_unmanaged/
#    Fig A1 — Ranked barplots: top 12 predictors, rows = model, cols = gap size
#    Fig A2 — Heatmap: rows = top predictors, cols = gap size, facet = model
#
#  Family C  →  01_ranked_barplots/ … 06_crosssite_summary/
#    Fig 1  — Ranked barplots per gap size (one per site × model × mgmt_var)
#    Fig 2  — Heatmap (one per site × mgmt_var, faceted by model)
#    Fig 3  — Rank-position distributions (RF only, median rank ± IQR ribbon)
#    Fig 4  — Management share of total importance (all models, faceted)
#    Fig 5  — Timelines: importance vs gap label index, RF, gap size M
#    Bonus  — Cross-site/run summary barplot (management vs top 8 others)
#
# PREDICTOR CLASSIFICATION
#   Each predictor is assigned a category via classify_predictor():
#     Management       — management variables and PI_{label} columns
#     Radiation        — PPFD, Rg, night
#     Meteorological   — Temp, VPD, RH, rain, rain_rolling_24
#     Phenology/time   — cyclical time encodings, season dummies
#     Other            — anything not matching the above
#   The category determines bar/line colour across all figures.
#
# PI COLLAPSE
#   The Phytomass Index produces one column per gap label (PI_L1, PI_M3, …).
#   These are collapsed to a single "PI" entry by averaging before any
#   plot is constructed, preventing the y-axis from being crowded with
#   identical per-label rows.
#
# INPUTS
#   results/{SITE}/{MANAGEMENT_OR_EFFECT}/{MODEL}/variable_importance/
#     rf_variable_importance_{SIZE}_all.rds  (aggregate, one per gap size)
#     rf_vi_{SIZE}_{LABEL}.rds               (per gap label, used by Figs 3, 5)
#
# OUTPUTS
#   graphs/variable_importance/
#     00_managed/, 00_unmanaged/  — Figs A1, A2
#     01_ranked_barplots/         — Fig 1
#     02_heatmaps/                — Fig 2
#     03_rank_distributions/      — Fig 3
#     04_mgmt_share/              — Fig 4
#     05_timelines/               — Fig 5
#     06_crosssite_summary/       — Bonus
#     session_info.txt
#
# SECTIONS
#    1  Packages & output directories
#    2  Configuration           — sites, models, management vars, colour palette
#    3  Generic VI loader       — load_vi_files() definition
#    4  Load Family A/B         — managed & unmanaged VI tables
#    5  Load Family C           — management effect VI tables
#    6  Shared helpers          — collapse PI labels, compute mean importance
#    7  Figures A1/B1           — ranked barplots for managed & unmanaged
#    8  Figures A2/B2           — heatmaps for managed & unmanaged
#    9  Figure 1                — ranked barplots for management effect
#   10  Figure 2                — heatmaps for management effect
#   11  Figure 3                — rank-position distributions (RF only)
#   12  Figure 4                — management share of importance
#   13  Figure 5                — timelines across gap labels (RF only)
#   14  Bonus                   — cross-site summary barplot
#   15  Session info
#
# HOW TO RUN
#   source("metrics/variable_importance_graphs.R")
#   The RF model runs must have completed first.
#
# =============================================================================


# =============================================================================
# SECTION 1 — Packages & output directories
# =============================================================================
# All output subdirectories are created here before any loop runs, so that
# save_fig() can write to any of them without further directory management.
# save_fig() saves both PDF and PNG versions of every figure.

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(forcats); library(stringr); library(glue)
})

OUT_DIR <- here::here("graphs", "variable_importance")

FIG_DIRS <- list(
  managed   = file.path(OUT_DIR, "00_managed"),
  unmanaged = file.path(OUT_DIR, "00_unmanaged"),
  fig1      = file.path(OUT_DIR, "01_ranked_barplots"),
  fig2      = file.path(OUT_DIR, "02_heatmaps"),
  fig3      = file.path(OUT_DIR, "03_rank_distributions"),
  fig4      = file.path(OUT_DIR, "04_mgmt_share"),
  fig5      = file.path(OUT_DIR, "05_timelines"),
  bonus     = file.path(OUT_DIR, "06_crosssite_summary")
)
invisible(lapply(FIG_DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

save_fig <- function(p, name, out_dir, w = 12, h = 7) {
  base <- file.path(out_dir, name)
  ggsave(paste0(base, ".pdf"), p, width = w, height = h)
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 300)
  message("  Saved: ", file.path(basename(out_dir), name))
}


# =============================================================================
# SECTION 2 — Configuration
# =============================================================================
# VI_PREFIX maps each model name to the prefix used in the VI filename.
# Currently only RF saves variable importance; extend this dictionary if
# MLP or XGBoost VI saving is added to their model scripts.
# classify_predictor() assigns categories based on exact column names;
# PI_{label} columns are caught by the str_starts("PI_") check.
# Add new predictors to the relevant case_when() arm as needed.

SITES     <- c("JC1", "JC2")
MODELS    <- c("RF", "MLP", "XGBoost")
GAP_SIZES <- c("VL", "L", "M", "S")

MGMT_VARS <- c(
  "Grazing_days_since",
  "Fertiliser_days_since",
  "N",
  "grass_height",
  "grass_biomass",
  "phytomass_index"
)

MGMT_LABELS <- c(
  "Grazing_days_since"    = "Days since grazing",
  "Fertiliser_days_since" = "Days since fertilisation",
  "N"                     = "N application",
  "grass_height"          = "Grass height",
  "grass_biomass"         = "Grass biomass",
  "phytomass_index"       = "Phytomass index (PI)"
)

# VI filename prefix per model — adjust if MLP/XGBoost use different naming
VI_PREFIX <- c(RF = "rf", MLP = "mlp", XGBoost = "xgb")

TOP_N <- 12

classify_predictor <- function(pred) {
  case_when(
    pred %in% c("Grazing_days_since", "Fertiliser_days_since",
                "N", "grass_height", "grass_biomass") |
      str_starts(pred, "PI_") | pred == "PI"            ~ "Management",
    pred %in% c("PPFD", "Rg", "night")                  ~ "Radiation",
    pred %in% c("Temp", "VPD", "RH", "rain",
                "rain_rolling_24")                       ~ "Meteorological",
    pred %in% c("hour_sin", "hour_cos",
                "doy_sin",  "doy_cos",
                "month_sin","month_cos",
                "Winter", "Spring", "Summer", "Autumn")  ~ "Phenology / time",
    TRUE                                                  ~ "Other"
  )
}

TYPE_COLOURS <- c(
  "Management"       = "#E07B39",
  "Radiation"        = "#F5C842",
  "Meteorological"   = "#4A90D9",
  "Phenology / time" = "#7BBF6A",
  "Other"            = "#AAAAAA"
)


# =============================================================================
# SECTION 3 — Generic VI loader — load_vi_files()
# =============================================================================
# Reads {prefix}_variable_importance_{SIZE}_all.rds for all four gap sizes
# from a given folder root.  Attaches site, model, gap_size, and any extra
# metadata columns passed via ... (e.g. management = "managed", mgmt_var = "N").
# Returns an empty data frame silently if the file does not exist, allowing
# the loading loops to skip missing combinations without error.

load_vi_files <- function(path_root, site, model, ...) {
  prefix  <- VI_PREFIX[[model]]
  extras  <- list(...)
  pmap_dfr(
    expand_grid(gap_size = GAP_SIZES),
    function(gap_size) {
      path <- file.path(
        path_root, "variable_importance",
        glue("{prefix}_variable_importance_{gap_size}_all.rds")
      )
      if (!file.exists(path)) return(NULL)
      df           <- readRDS(path)
      df$site      <- site
      df$model     <- model
      df$gap_size  <- gap_size
      for (nm in names(extras)) df[[nm]] <- extras[[nm]]
      df
    }
  )
}


# =============================================================================
# SECTION 4 — Load Family A/B: managed and unmanaged VI tables
# =============================================================================
# Iterates over management conditions × sites × models and stacks all
# available VI files into vi_base.  The management column distinguishes
# managed (full predictor set) from unmanaged (BASE only) runs.
# Predictor type and is_mgmt flag are added here for use in all downstream
# figure functions.

message("Loading managed & unmanaged VI files ...")

vi_base <- map_dfr(c("managed", "unmanaged"), function(mgmt) {
  map_dfr(SITES, function(s) {
    map_dfr(MODELS, function(m) {
      root <- here::here("results", s, mgmt, m)
      load_vi_files(root, site = s, model = m, management = mgmt)
    })
  })
}) %>%
  mutate(
    pred_type  = classify_predictor(predictor),
    is_mgmt    = pred_type == "Management",
    gap_size   = factor(gap_size, levels = GAP_SIZES),
    management = factor(management, levels = c("unmanaged", "managed"))
  )

message("  Rows loaded (managed + unmanaged): ", nrow(vi_base))


# =============================================================================
# SECTION 5 — Load Family C: management effect VI tables
# =============================================================================
# Iterates over sites × management variables × models.  Each management
# variable has its own subfolder under management_effect/RF/ (or MLP/, XGBoost/).
# The mgmt_var column in vi_mgmt_effect identifies which ablation run each
# row came from.

message("Loading management effect VI files ...")

vi_mgmt_effect <- pmap_dfr(
  expand_grid(site = SITES, mgmt_var = MGMT_VARS, model = MODELS),
  function(site, mgmt_var, model) {
    root <- here::here("results", site, "management_effect", model, mgmt_var)
    df   <- load_vi_files(root, site = site, model = model, mgmt_var = mgmt_var)
    df
  }
) %>%
  mutate(
    pred_type = classify_predictor(predictor),
    is_mgmt   = pred_type == "Management",
    gap_size  = factor(gap_size, levels = GAP_SIZES)
  )

message("  Rows loaded (management effect): ", nrow(vi_mgmt_effect))


# =============================================================================
# SECTION 6 — Shared helper: collapse PI labels and compute mean importance
# =============================================================================
# collapse_and_mean() performs two operations:
#   1. Renames all PI_{label} predictors (PI_L1, PI_M3, …) to "PI" so that
#      all per-gap PI columns are treated as one predictor in plots.
#      Averaging across the renamed rows gives the mean PI importance.
#   2. Computes mean importance_rel across gap labels within each group.
# The group_vars argument specifies which columns define the grouping
# (e.g. site × management × model × gap_size for the base runs).

collapse_and_mean <- function(df, group_vars) {
  df %>%
    mutate(predictor = if_else(str_starts(predictor, "PI_"), "PI", predictor)) %>%
    group_by(across(all_of(c(group_vars, "predictor", "pred_type", "is_mgmt")))) %>%
    summarise(mean_imp = mean(importance_rel, na.rm = TRUE), .groups = "drop")
}

# Pre-compute mean importance for both families
vi_base_mean   <- collapse_and_mean(vi_base,
                                    c("site", "management", "model", "gap_size"))
vi_effect_mean <- collapse_and_mean(vi_mgmt_effect,
                                    c("site", "model", "mgmt_var", "gap_size"))

# Management predictor names after PI collapse (for heatmap filter)
mgmt_preds <- vi_mgmt_effect %>%
  filter(is_mgmt) %>%
  mutate(predictor = if_else(str_starts(predictor, "PI_"), "PI", predictor)) %>%
  distinct(predictor) %>%
  pull(predictor)


# =============================================================================
# SECTION 7 — Figures A1/B1: Ranked barplots — managed and unmanaged
# =============================================================================
# Produces one horizontal bar chart per site × management condition.
# Layout: rows = models (RF, MLP, XGBoost), columns = gap sizes (VL, L, M, S).
# Each panel shows the top TOP_N predictors ranked by mean relative importance,
# with bars coloured by predictor type.  This figure answers the question:
# "Which predictors dominate the model when (no) management data is available?"

message("\nFig A1/B1 — Ranked barplots: managed & unmanaged ...")

for (mgmt in c("managed", "unmanaged")) {
  for (s in SITES) {
    dat <- vi_base_mean %>%
      filter(site == s, management == mgmt) %>%
      group_by(model, gap_size) %>%
      slice_max(mean_imp, n = TOP_N) %>%
      ungroup() %>%
      mutate(predictor = fct_reorder2(predictor, gap_size, mean_imp))
    
    if (nrow(dat) == 0) next
    
    p <- ggplot(dat, aes(x = mean_imp, y = predictor, fill = pred_type)) +
      geom_col(width = 0.7) +
      facet_grid(model ~ gap_size, scales = "free_y",
                 labeller = labeller(gap_size = c(
                   VL = "Very long", L = "Long", M = "Medium", S = "Short"
                 ))) +
      scale_fill_manual(values = TYPE_COLOURS, name = "Predictor type") +
      scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      labs(
        title    = glue("Top {TOP_N} predictors — {s} | {str_to_title(mgmt)}"),
        subtitle = "Mean relative importance across all artificial gap labels",
        x = "Mean relative importance (%)", y = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom",
            strip.background = element_rect(fill = "grey92"),
            panel.grid.minor = element_blank())
    
    save_fig(p, glue("barplot_{s}"), out_dir = FIG_DIRS[[mgmt]], w = 14, h = 10)
  }
}


# =============================================================================
# SECTION 8 — Figures A2/B2: Heatmaps — managed and unmanaged
# =============================================================================
# For each site × management condition, displays the importance of the top 10
# predictors (by overall average across models and gap sizes) as a heatmap.
# Rows = predictors, columns = gap sizes, faceted by model.
# Cell colour encodes mean relative importance using the inferno scale (dark =
# low, yellow = high).  Management variable rows are marked with ★ and bolded.

message("\nFig A2/B2 — Heatmaps: managed & unmanaged ...")

for (mgmt in c("managed", "unmanaged")) {
  for (s in SITES) {
    dat_all <- vi_base_mean %>% filter(site == s, management == mgmt)
    if (nrow(dat_all) == 0) next
    
    top_preds <- dat_all %>%
      group_by(predictor) %>%
      summarise(overall = mean(mean_imp, na.rm = TRUE)) %>%
      slice_max(overall, n = 10) %>%
      pull(predictor)
    
    dat <- dat_all %>%
      filter(predictor %in% top_preds) %>%
      mutate(
        pred_label = if_else(is_mgmt, paste0("★ ", predictor), predictor),
        pred_label = fct_reorder(pred_label, mean_imp)
      )
    
    p <- ggplot(dat, aes(x = gap_size, y = pred_label, fill = mean_imp)) +
      geom_tile(colour = "white", linewidth = 0.4) +
      geom_text(aes(label = scales::percent(mean_imp, accuracy = 0.1)),
                size = 2.8, colour = "white") +
      facet_wrap(~model, nrow = 1) +
      scale_fill_viridis_c(option = "inferno", direction = -1,
                           labels = scales::percent_format(accuracy = 0.1),
                           name   = "Mean relative\nimportance (%)") +
      scale_x_discrete(labels = c(VL = "Very long", L = "Long",
                                  M = "Medium",    S = "Short")) +
      labs(
        title    = glue("Variable importance heatmap — {s} | {str_to_title(mgmt)}"),
        subtitle = "Top predictors averaged across all gap labels  |  ★ = management",
        x = "Gap size", y = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "right", panel.grid = element_blank(),
            strip.background = element_rect(fill = "grey92"))
    
    save_fig(p, glue("heatmap_{s}"), out_dir = FIG_DIRS[[mgmt]], w = 13, h = 7)
  }
}


# =============================================================================
# SECTION 9 — Figure 1: Ranked barplots — management effect
# =============================================================================
# Produces one figure per site × model × management variable.
# Each figure shows the top TOP_N predictors in a 2×2 facet grid of gap sizes.
# The management variable's bar reveals where it ranks relative to the
# meteorological predictors.  If it falls outside the top TOP_N, it does not
# appear — this is a scientifically meaningful result (low importance).

message("\nFig 1 — Ranked barplots: management effect ...")

for (s in SITES) {
  for (mv in MGMT_VARS) {
    for (m in MODELS) {
      dat <- vi_effect_mean %>%
        filter(site == s, mgmt_var == mv, model == m) %>%
        group_by(gap_size) %>%
        slice_max(mean_imp, n = TOP_N) %>%
        ungroup() %>%
        mutate(predictor = fct_reorder2(predictor, gap_size, mean_imp))
      
      if (nrow(dat) == 0) next
      
      p <- ggplot(dat, aes(x = mean_imp, y = predictor, fill = pred_type)) +
        geom_col(width = 0.7) +
        facet_wrap(~gap_size, scales = "free_y", ncol = 2,
                   labeller = labeller(gap_size = c(
                     VL = "Very long gap (VL)", L = "Long gap (L)",
                     M  = "Medium gap (M)",     S = "Short gap (S)"
                   ))) +
        scale_fill_manual(values = TYPE_COLOURS, name = "Predictor type") +
        scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
        labs(
          title    = glue("Top {TOP_N} predictors — {s} | {m} | {MGMT_LABELS[mv]}"),
          subtitle = "Mean importance_rel across all artificial gap labels",
          x = "Mean relative importance (%)", y = NULL
        ) +
        theme_bw(base_size = 11) +
        theme(legend.position = "bottom",
              strip.background = element_rect(fill = "grey92"),
              panel.grid.minor = element_blank())
      
      save_fig(p, glue("fig1_barplot_{s}_{m}_{mv}"),
               out_dir = FIG_DIRS$fig1, w = 11, h = 9)
    }
  }
}


# =============================================================================
# SECTION 10 — Figure 2: Heatmaps — management effect
# =============================================================================
# Rows = management predictors + top 6 non-management predictors (by overall
# average importance across all ablation runs).  Columns = gap sizes.
# Faceted by model.  One figure per site × management variable.
# Management predictor rows are marked with ★ and bolded on the y-axis.
# This figure is designed to show, for each model simultaneously, how the
# added management variable ranks and whether it displaces any meteorological
# predictors.

message("\nFig 2 — Heatmaps: management effect ...")

top_nonmgmt_effect <- vi_effect_mean %>%
  filter(!is_mgmt) %>%
  group_by(predictor) %>%
  summarise(overall = mean(mean_imp, na.rm = TRUE)) %>%
  slice_max(overall, n = 6) %>%
  pull(predictor)

heatmap_preds <- c(mgmt_preds, top_nonmgmt_effect)

for (s in SITES) {
  for (mv in MGMT_VARS) {
    dat <- vi_effect_mean %>%
      filter(site == s, mgmt_var == mv,
             predictor %in% heatmap_preds) %>%
      mutate(
        pred_label = if_else(is_mgmt, paste0("★ ", predictor), predictor),
        pred_label = fct_reorder(pred_label, mean_imp)
      )
    
    if (nrow(dat) == 0) next
    
    p <- ggplot(dat, aes(x = gap_size, y = pred_label, fill = mean_imp)) +
      geom_tile(colour = "white", linewidth = 0.4) +
      geom_text(aes(label = scales::percent(mean_imp, accuracy = 0.1)),
                size = 2.8, colour = "white") +
      facet_wrap(~model, nrow = 1) +
      scale_fill_viridis_c(option = "inferno", direction = -1,
                           labels = scales::percent_format(accuracy = 0.1),
                           name   = "Mean relative\nimportance (%)") +
      scale_x_discrete(labels = c(VL = "Very long", L = "Long",
                                  M = "Medium",    S = "Short")) +
      labs(
        title    = glue("Variable importance heatmap — {s} | {MGMT_LABELS[mv]}"),
        subtitle = "★ = management predictors  |  facets = model",
        x = "Gap size", y = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.y      = element_text(
          face = if_else(str_starts(levels(dat$pred_label), "★"), "bold", "plain")
        ),
        legend.position  = "right",
        panel.grid       = element_blank(),
        strip.background = element_rect(fill = "grey92")
      )
    
    save_fig(p, glue("fig2_heatmap_{s}_{mv}"),
             out_dir = FIG_DIRS$fig2, w = 12, h = 7)
  }
}


# =============================================================================
# SECTION 11 — Figure 3: Rank-position distributions (RF only)
# =============================================================================
# For each gap label, each predictor is ranked by its absolute importance
# (rank 1 = most important).  The median rank and IQR across all gap labels
# are then plotted vs gap size on a reversed y-axis (rank 1 at top).
# An IQR ribbon shows the spread; management predictor lines are solid,
# meteorological lines are dashed.  A narrow ribbon indicates that the
# predictor's rank is stable across gap labels; a wide ribbon indicates
# that its importance is highly variable.

message("\nFig 3 — Rank-position distributions (RF) ...")

vi_ranked_effect <- vi_mgmt_effect %>%
  filter(model == "RF") %>%
  group_by(site, mgmt_var, gap_size, gap_label) %>%
  mutate(rank_pos = rank(-importance, ties.method = "min")) %>%
  ungroup()

vi_rank_summary <- vi_ranked_effect %>%
  group_by(site, mgmt_var, gap_size, predictor, pred_type, is_mgmt) %>%
  summarise(
    median_rank = median(rank_pos, na.rm = TRUE),
    q25         = quantile(rank_pos, 0.25, na.rm = TRUE),
    q75         = quantile(rank_pos, 0.75, na.rm = TRUE),
    .groups     = "drop"
  )

top_nonmgmt_rank <- vi_rank_summary %>%
  filter(!is_mgmt) %>%
  group_by(predictor) %>%
  summarise(med = median(median_rank)) %>%
  slice_min(med, n = 5) %>%
  pull(predictor)

for (s in SITES) {
  for (mv in MGMT_VARS) {
    dat <- vi_rank_summary %>%
      filter(site == s, mgmt_var == mv,
             predictor %in% c(mgmt_preds, top_nonmgmt_rank)) %>%
      mutate(predictor = fct_reorder(predictor, median_rank))
    
    if (nrow(dat) == 0) next
    
    p <- ggplot(dat, aes(x = gap_size, y = median_rank,
                         colour = pred_type, group = predictor)) +
      geom_line(aes(linetype = is_mgmt), linewidth = 0.9) +
      geom_ribbon(aes(ymin = q25, ymax = q75, fill = pred_type),
                  alpha = 0.15, colour = NA) +
      geom_point(size = 2.5) +
      scale_colour_manual(values = TYPE_COLOURS, name = "Predictor type") +
      scale_fill_manual(values = TYPE_COLOURS, guide = "none") +
      scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dashed"),
                            labels = c("TRUE" = "Management", "FALSE" = "Other"),
                            name = NULL) +
      scale_y_reverse(name = "Median importance rank (1 = top)") +
      scale_x_discrete(labels = c(VL = "Very long", L = "Long",
                                  M = "Medium",    S = "Short")) +
      labs(
        title    = glue("Rank-position — {s} | RF | {MGMT_LABELS[mv]}"),
        subtitle = "Ribbon = IQR across gap labels.  Lower = more important.",
        x = "Gap size"
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank())
    
    save_fig(p, glue("fig3_rank_{s}_{mv}"),
             out_dir = FIG_DIRS$fig3, w = 10, h = 6)
  }
}


# =============================================================================
# SECTION 12 — Figure 4: Management share of total importance
# =============================================================================
# For each gap label, the management share is sum(importance_rel[is_mgmt]).
# Boxplots across gap labels show the distribution of this share per management
# variable × model × gap size.  A high and stable management share across gap
# sizes and models indicates that the management signal is consistently captured
# by the model, not driven by a few atypical gaps.

message("\nFig 4 — Management share of total importance ...")

vi_mgmt_share <- vi_mgmt_effect %>%
  group_by(site, model, mgmt_var, gap_size, gap_label) %>%
  summarise(
    mgmt_share = sum(importance_rel[is_mgmt], na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  mutate(gap_size = factor(gap_size, levels = GAP_SIZES))

for (s in SITES) {
  dat <- vi_mgmt_share %>% filter(site == s)
  if (nrow(dat) == 0) next
  
  p <- ggplot(dat, aes(x = gap_size, y = mgmt_share, fill = mgmt_var)) +
    geom_boxplot(alpha = 0.75, outlier.size = 1.2, width = 0.55,
                 position = position_dodge(width = 0.7)) +
    facet_wrap(~model, nrow = 1) +
    scale_fill_brewer(palette = "Set2",
                      labels  = MGMT_LABELS[unique(dat$mgmt_var)],
                      name    = "Management variable") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                       name   = "Management share of total importance (%)") +
    scale_x_discrete(labels = c(VL = "Very long", L = "Long",
                                M = "Medium",    S = "Short")) +
    labs(
      title    = glue("Management share of importance — {s}"),
      subtitle = "Distribution across gap labels  |  facets = model",
      x = "Gap size"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92"))
  
  save_fig(p, glue("fig4_mgmt_share_{s}"),
           out_dir = FIG_DIRS$fig4, w = 13, h = 6)
}


# =============================================================================
# SECTION 13 — Figure 5: Importance timelines (RF only, gap size M)
# =============================================================================
# Plots the relative importance of the management predictor(s) against the
# sequential gap label index (a proxy for time within the two-year record).
# A LOESS smoother with 95% confidence band highlights temporal trends.
# Seasonal patterns in importance (e.g. higher in summer for grass biomass,
# higher in spring for grazing days since) would support the ecological
# interpretation of the predictor's contribution.
# Gap size M (medium, ~7 days) is used by default as it provides enough gap
# labels for a meaningful temporal trend while avoiding the very sparse
# sampling of the VL category.  Change TIMELINE_SIZE to use a different size.

message("\nFig 5 — Management variable timelines (RF, gap size M) ...")

TIMELINE_SIZE <- "M"

vi_timeline <- vi_mgmt_effect %>%
  filter(model == "RF", gap_size == TIMELINE_SIZE, is_mgmt) %>%
  mutate(gap_index = as.integer(str_extract(gap_label, "\\d+")))

for (s in SITES) {
  for (mv in MGMT_VARS) {
    dat <- vi_timeline %>% filter(site == s, mgmt_var == mv)
    if (nrow(dat) == 0) next
    
    # Each PI_{label} is a distinct line — collapse to "PI" for cleaner legend
    dat <- dat %>%
      mutate(predictor = if_else(str_starts(predictor, "PI_"), "PI", predictor)) %>%
      group_by(site, mgmt_var, gap_index, predictor, pred_type, is_mgmt) %>%
      summarise(importance_rel = mean(importance_rel, na.rm = TRUE), .groups = "drop")
    
    p <- ggplot(dat, aes(x = gap_index, y = importance_rel,
                         colour = predictor, group = predictor)) +
      geom_line(linewidth = 0.8, alpha = 0.8) +
      geom_point(size = 2, alpha = 0.8) +
      geom_smooth(se = TRUE, method = "loess", span = 0.4,
                  linewidth = 0.5, alpha = 0.12) +
      scale_colour_manual(
        values = setNames(
          rep(TYPE_COLOURS[["Management"]], length(unique(dat$predictor))),
          unique(dat$predictor)
        ),
        name = "Predictor"
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      labs(
        title    = glue("Importance over time — {s} | RF | {MGMT_LABELS[mv]} | gap size: {TIMELINE_SIZE}"),
        subtitle = "Each point = one artificial gap label (ordered in time).  Shaded band = LOESS ± SE.",
        x = glue("Gap label index ({TIMELINE_SIZE}1 → {TIMELINE_SIZE}N)"),
        y = "Relative importance (%)"
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank())
    
    save_fig(p, glue("fig5_timeline_{s}_{mv}"),
             out_dir = FIG_DIRS$fig5, w = 11, h = 5)
  }
}


# =============================================================================
# SECTION 14 — Bonus: Cross-site/model/run summary barplot
# =============================================================================
# Aggregates mean importance across all sites, models, gap sizes, and
# management ablation runs to produce a single "headline" barplot.
# Shows management variables alongside the top 8 non-management predictors,
# with error bars of ±1 SE across the aggregation groups.  This figure
# provides the most compact summary of the overall importance landscape and
# is suitable as a main-text figure.

message("\nBonus — Cross-site summary barplot ...")

vi_cross <- vi_effect_mean %>%
  group_by(predictor, pred_type, is_mgmt) %>%
  summarise(
    mean_imp_all = mean(mean_imp, na.rm = TRUE),
    se_imp       = sd(mean_imp,   na.rm = TRUE) / sqrt(n()),
    .groups      = "drop"
  ) %>%
  {
    top8_non <- filter(., !is_mgmt) %>%
      slice_max(mean_imp_all, n = 8) %>% pull(predictor)
    filter(., is_mgmt | predictor %in% top8_non)
  } %>%
  mutate(predictor = fct_reorder(predictor, mean_imp_all))

p_cross <- ggplot(vi_cross,
                  aes(x = mean_imp_all, y = predictor, fill = pred_type)) +
  geom_col(width = 0.7) +
  geom_errorbarh(aes(xmin = mean_imp_all - se_imp,
                     xmax = mean_imp_all + se_imp),
                 height = 0.3, linewidth = 0.5, colour = "grey40") +
  scale_fill_manual(values = TYPE_COLOURS, name = "Predictor type") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    title    = "Mean relative importance — all sites, models, gap sizes, management runs",
    subtitle = "Error bars = ± 1 SE.  Management predictors highlighted in orange.",
    x = "Mean relative importance (%)", y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

save_fig(p_cross, "fig_bonus_crosssite_summary",
         out_dir = FIG_DIRS$bonus, w = 10, h = 7)


# =============================================================================
# SECTION 15 — Session info
# =============================================================================
# Writes R session information (package versions, platform) to a text file
# in the output directory for reproducibility documentation.

message("\n", strrep("=", 60))
message("All figures saved under: ", OUT_DIR)
message(strrep("=", 60))
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "session_info.txt"))