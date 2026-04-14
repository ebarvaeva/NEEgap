# =============================================================================
# plot_vi_summary.R — RF Variable Importance Summary Boxplots
# =============================================================================
#
# PURPOSE
#   Produces publication-ready figures and CSV summaries of RF variable
#   importance across sites (JC1, JC2), gap sizes (S, M, L, VL), and
#   predictor sets (managed baseline and feature-addition runs).
#
#   Variable importance is read directly from pre-computed RDS files —
#   no model re-fitting is performed here.
#
# =============================================================================
#
# INPUT FILES
#   Plot 1 — Managed predictor set:
#     results/{SITE}/managed/RF/variable_importance/
#       rf_variable_importance_S_all.rds
#       rf_variable_importance_M_all.rds
#       rf_variable_importance_L_all.rds
#       rf_variable_importance_VL_all.rds
#
#   Plot 2 — Feature-addition runs (BASE + one management variable at a time):
#     results/{SITE}/management_effect/RF/{MGMT_VAR}/variable_importance/
#       rf_variable_importance_S_all.rds
#       rf_variable_importance_M_all.rds
#       rf_variable_importance_L_all.rds
#       rf_variable_importance_VL_all.rds
#
#   Each RDS file contains a data frame with at least two columns:
#     predictor      — name of the input variable
#     importance_rel — relative importance in [0, 1], where
#                      sum across all predictors = 1 for each gap Si
#
#   PI columns are stored as PI_{gap_label} (e.g. PI_S1, PI_M3) in both
#   the managed run and the phytomass_index feature-addition run. Both are
#   collapsed to the unified label "PI" before plotting.
#
# =============================================================================
#
# OUTPUT FILES
#   graphs/VI/plot1_managed_RF.{png,pdf}
#   graphs/VI/plot2_{MGMT_VAR}_RF.{png,pdf}   — one per management variable
#     e.g. plot2_Grazing_days_since_RF.png
#          plot2_Fertiliser_days_since_RF.png
#          plot2_N_RF.png
#          plot2_grass_height_RF.png
#          plot2_grass_biomass_RF.png
#          plot2_phytomass_index_RF.png
#
#   graphs/VI/csv/vi_summary_managed.csv
#   graphs/VI/csv/vi_summary_{MGMT_VAR}.csv   — one per management variable
#
# =============================================================================
#
# PREDICTOR FILTERING
#   Only predictors satisfying at least one of the following are plotted:
#     (a) median importance_rel > IMP_THRESHOLD (default 0) across all
#         gaps, gap sizes, and sites within the run
#     (b) the predictor is one of the 6 management variables (always shown
#         regardless of importance, so their contribution is always visible)
#
# =============================================================================
#
# PREDICTOR LABELS
#   y-axis labels follow the Symbol/name column of Table tab:variables.
#   Code names are mapped to display symbols before plotting.
#   Subscripts use R plotmath notation (parsed by scale_x_discrete):
#     T[air]  →  T subscript air
#     R[g]    →  R subscript g
#     u['*']  →  u subscript *
#     dsg     →  dsg (no subscript)
#     dsf     →  dsf (no subscript)
#     h[s]    →  h subscript s
#   Any code name not in the lookup is shown as-is (fallback).
#
# =============================================================================
#
# PLOT 1 — Managed Predictor Set
# --------------------------------
# Shows the full distribution of relative importance for every retained
# predictor in the managed predictor set, pooled across all artificial gaps
# and gap sizes (S, M, L, VL).
#
# Layout:
#   y-axis  : predictor symbol (Table tab:variables), ordered top → bottom
#              by descending median importance across all gap sizes and both sites
#   x-axis  : relative importance [0, 1]
#   Boxes   : 2 boxplots per predictor, one per site (JC1 = white, JC2 = black),
#             pooling all gap sizes (S / M / L / VL)
#
# Pipeline:
#   Step 1 — For gap size G ∈ {S, M, L, VL}, read the RDS file containing
#            K_G individual gap importance values (one row per gap Si):
#              importance_rel(predictor p, gap Si, size G)
#
#   Step 2 — Collapse all PI_{Si} rows to the unified label "PI":
#              if predictor starts with "PI_"     →  rename to "PI"
#              if predictor == "phytomass_index"  →  rename to "PI"
#            Each PI_{Si} row becomes a separate "PI" data point, contributing
#            individually to the PI boxplot IQR — this is correct.
#
#   Step 3 — Filter predictors to keep:
#              median{ importance_rel(p) } > IMP_THRESHOLD  OR
#              p ∈ {management variable labels}
#
#   Step 4 — For each retained predictor p and site, the boxplot pools all
#            gaps across all gap sizes and shows:
#              Median   = median{ importance_rel(p, Si, G) : all i, all G }
#              IQR      = Q75 - Q25 of the same distribution
#              Whiskers = 1.5 × IQR beyond Q25 / Q75
#              Outliers = individual gap values beyond the whiskers
#
#   Step 5 — Predictor order on y-axis:
#              rank(p) = median{ importance_rel(p, Si, G, site) :
#                                all i, all G ∈ {S,M,L,VL}, both sites }
#            The most important predictor appears at the top.
#
#   Step 6 — y-axis labels mapped to table symbols via PREDICTOR_LABELS.
#
# =============================================================================
#
# PLOT 2 — Feature-Addition Runs (one plot per management variable)
# -----------------------------------------------------------------
# For each management variable m, produces one plot showing the full predictor
# importance distribution from the BASE+m feature-addition run, pooled across
# all gap sizes. Structure is identical to Plot 1.
#
# Layout (per plot, one per MGMT_VAR):
#   y-axis  : retained predictor symbols, ordered top → bottom by descending
#              median importance within this run
#   x-axis  : relative importance [0, 1]
#   Boxes   : 2 boxplots per predictor (JC1 = white, JC2 = black)
#   Subtitle: "BASE + {variable label}" identifies which run is shown
#
# Pipeline:
#   Step 1 — Read RDS files and collapse PI names (same as Plot 1).
#   Step 2 — Filter: keep median > IMP_THRESHOLD OR management variable.
#   Step 3 — Boxplot IQR pooled from all S1…Sn gaps across all gap sizes.
#   Step 4 — Predictor order by descending median within this run.
#   Step 5 — y-axis labels mapped to table symbols via PREDICTOR_LABELS.
#
# =============================================================================
#
# PI NAMING CONVENTION
#   In both the managed run and the phytomass_index feature-addition run, the
#   Phytomass Index is recomputed per gap and stored as separate columns
#   PI_S1, PI_S2, PI_M1, ... (one per gap label). After tidy_vi_raw(),
#   all of these become "PI" and each contributes one data point to the
#   "PI" boxplot — one observation per gap, which is correct.
#
# =============================================================================
#
# HOW TO RUN
#   Ensure results/ and graphs/ directories exist at the project root, then:
#   source("scripts/plot_vi_summary.R")
#
# DEPENDENCIES
#   here, dplyr, tidyr, purrr, ggplot2, forcats, stringr, glue, scales
# =============================================================================
# =============================================================================
# SECTION 1 — Settings
# =============================================================================

RESULTS_ROOT <- here::here("results")
OUT_DIR      <- here::here("graphs", "VI")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SITES     <- c("JC1", "JC2")
GAP_SIZES <- c("S", "M", "L", "VL")

# Importance threshold — predictors with median importance_rel <= this value
# are dropped UNLESS they are one of the 6 management variables
IMP_THRESHOLD <- 0

# Management variables for Plot 2 — order controls plot production sequence.
# These are ALWAYS shown in plots regardless of importance threshold.
MGMT_VARS <- c(
  "Grazing_days_since",
  "Fertiliser_days_since",
  "N",
  "grass_height",
  "grass_biomass",
  "phytomass_index"
)

# Labels checked against always-show list (use unified "PI" after tidy_vi_raw)
MGMT_PLOT_LABELS <- c(
  "Grazing_days_since",
  "Fertiliser_days_since",
  "N",
  "grass_height",
  "grass_biomass",
  "PI"
)

# Short display labels used in plot subtitles
MGMT_COL_LABELS <- c(
  Grazing_days_since    = "Grazing",
  Fertiliser_days_since = "Fertiliser",
  N                     = "N",
  grass_height          = "Height",
  grass_biomass         = "Biomass",
  phytomass_index       = "PI"
)

# Predictor display labels — maps code names to Symbol/name from Table tab:variables.
# Uses R plotmath notation; parsed by scale_x_discrete in make_vi_boxplot().
# Subscripts: X[y] = X subscript y.
# Any code name not listed here falls back to the raw code name.
PREDICTOR_LABELS <- c(
  # Meteorological and radiation drivers
  PPFD              = "PPFD",
  night             = "night",
  ustar             = "u['*']",          # u*
  Temp              = "T[air]",          # T_air
  RH                = "RH",
  rain              = "P",               # P
  rain_rolling_24   = "P[24]",           # P_24
  Rg                = "R[g]",            # R_g
  VPD               = "VPD",
  # Temporal cyclic encodings
  hour_sin          = "hour_sin",
  hour_cos          = "hour_cos",
  doy_sin           = "doy_sin",
  doy_cos           = "doy_cos",
  month_sin         = "month_sin",
  month_cos         = "month_cos",
  # Seasonal indicators
  Winter            = "Winter",
  Spring            = "Spring",
  Summer            = "Summer",
  Autumn            = "Autumn",
  # Management variables
  # Management variables
  Grazing_days_since    = "bold(dsg)",
  Fertiliser_days_since = "bold(dsf)",
  grass_height          = "bold(h[s])",
  grass_biomass         = "bold(B)",
  N                     = "bold(N)",
  PI                    = "bold(PI)"
)

# Colours per gap size — blue family for short gaps, red family for long
GAP_COLOURS <- c(
  S  = "#4575B4",   # deep blue
  M  = "#74ADD1",   # light blue
  L  = "#FDB863",   # light orange
  VL = "#D73027"    # deep red
)

GAP_LABELS <- c(S = "S", M = "M", L = "L", VL = "VL")


# =============================================================================
# SECTION 2 — Packages
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(dplyr); library(tidyr); library(purrr)
  library(ggplot2); library(forcats); library(stringr); library(glue)
  library(scales)
})

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p,
         width = w, height = h)
  message("Saved: ", name)
}


# =============================================================================
# SECTION 3 — VI loader, tidier, and filter
# =============================================================================

# load_rf_vi_raw()
#   Reads rf_variable_importance_{G}_all.rds for all G ∈ {S, M, L, VL}
#   from a given folder. Returns all individual gap rows — NOT averaged.
#   Each row = one gap Si's importance value for one predictor.

load_rf_vi_raw <- function(vi_root) {
  map_dfr(GAP_SIZES, function(gs) {
    path <- file.path(vi_root, "variable_importance",
                      glue("rf_variable_importance_{gs}_all.rds"))
    if (!file.exists(path)) {
      message("  Not found (skip): ", path)
      return(NULL)
    }
    readRDS(path) %>%
      mutate(gap_size = gs)
  })
}

# tidy_vi_raw()
#   Collapses all PI-related predictor names to the unified label "PI":
#     PI_S1, PI_S2, PI_M1, ...  →  "PI"   (per-gap columns in managed run)
#     phytomass_index            →  "PI"   (ablation run naming)
#     PI (already plain)         →  "PI"   (no change)
#   All other predictors are unchanged.

tidy_vi_raw <- function(df) {
  df %>%
    mutate(predictor = case_when(
      str_starts(predictor, "PI_")   ~ "PI",
      predictor == "phytomass_index" ~ "PI",
      TRUE                           ~ predictor
    )) %>%
    mutate(gap_size = factor(gap_size, levels = GAP_SIZES))
}

# filter_predictors()
#   Keeps predictor p if:
#     (a) median importance_rel across all rows > IMP_THRESHOLD, OR
#     (b) p is one of the 6 management variable labels (always shown)

filter_predictors <- function(df) {
  med_imp <- df %>%
    group_by(predictor) %>%
    summarise(med = median(importance_rel, na.rm = TRUE), .groups = "drop")
  
  keep <- med_imp %>%
    filter(med > IMP_THRESHOLD | predictor %in% MGMT_PLOT_LABELS) %>%
    pull(predictor)
  
  dropped <- setdiff(med_imp$predictor, keep)
  if (length(dropped))
    message("  Dropped (median <= ", IMP_THRESHOLD, "): ",
            paste(dropped, collapse = ", "))
  
  df %>% filter(predictor %in% keep)
}

# make_vi_boxplot()
#   Shared plot builder used by both Plot 1 and Plot 2.
#   x-axis labels are mapped from code names to table symbols via
#   PREDICTOR_LABELS, parsed as plotmath expressions for subscripts.
#   Any code name not in PREDICTOR_LABELS is shown as-is (fallback).

make_vi_boxplot <- function(df, pred_order, subtitle = NULL) {
  
  df <- df %>%
    mutate(
      predictor = factor(predictor, levels = pred_order),
      site      = factor(site,      levels = SITES),
      gap_size  = factor(gap_size,  levels = GAP_SIZES)
    )
  
  x_labels <- function(x) {
    labs <- PREDICTOR_LABELS[x]
    labs[is.na(labs)] <- x[is.na(labs)]
    parse(text = labs)
  }
  
  SITE_COLOURS <- c(JC1 = "blue", JC2 = "red")
  
  p <- ggplot(df,
              aes(x    = predictor,
                  y    = importance_rel,
                  fill = site)) +             # ← fill encodes site
    geom_boxplot(
      outlier.size  = 0.9,
      outlier.alpha = 0.9,
      linewidth     = 0.3,
      position      = position_dodge(width = 0.9),   # ← side-by-side
      width         = 0.9
    ) +
    coord_flip() +
    scale_fill_manual(values = SITE_COLOURS, name = "Site") +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(
      name   = "Relative importance",
      breaks = pretty_breaks(n = 6),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    labs(x = NULL, subtitle = subtitle) +
    theme_bw(base_size = 35) +
    theme(
      legend.position    = "bottom",
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background   = element_rect(fill = "grey92"),
      strip.text         = element_text(size = 32, face = "bold"),
      axis.text.x        = element_text(angle = 32, hjust = 1, size = 32),
      axis.text.y        = element_text(size = 32),
      axis.title.y       = element_text(size = 32),
      legend.text        = element_text(size = 32),
      legend.title       = element_text(size = 32),
      plot.subtitle      = element_text(size = 32, face = "bold")
    )
  p
}

# =============================================================================
# SECTION 4 — PLOT 1: Managed predictor set
# =============================================================================

message("\n[Plot 1] Loading managed RF VI files ...")

vi_managed_raw <- map_dfr(SITES, function(s) {
  vi_root <- file.path(RESULTS_ROOT, s, "managed", "RF")
  raw     <- load_rf_vi_raw(vi_root)
  if (is.null(raw) || nrow(raw) == 0) {
    message("  Not found: ", vi_root); return(NULL)
  }
  tidy_vi_raw(raw) %>% mutate(site = s)
})

if (nrow(vi_managed_raw) == 0) stop("No managed VI data found.")

vi_managed_filtered <- filter_predictors(vi_managed_raw)

predictor_order_managed <- vi_managed_filtered %>%
  group_by(predictor) %>%
  summarise(med = median(importance_rel, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(predictor)

p1 <- make_vi_boxplot(vi_managed_filtered, predictor_order_managed,
                      subtitle = "")

n_preds <- length(predictor_order_managed)

save_fig(p1, "plot1_managed_RF", w = 30, h = 17)

# =============================================================================
# SECTION 5 — PLOT 2: Management effect ablation runs
# =============================================================================

message("\n[Plot 2] Loading management effect RF VI files ...")

for (mv in MGMT_VARS) {
  
  message("  Processing: ", mv)
  
  vi_mv_raw <- map_dfr(SITES, function(s) {
    vi_root <- file.path(RESULTS_ROOT, s, "management_effect", "RF", mv)
    raw     <- load_rf_vi_raw(vi_root)
    if (is.null(raw) || nrow(raw) == 0) {
      message("    Not found: ", vi_root); return(NULL)
    }
    tidy_vi_raw(raw) %>% mutate(site = s)
  })
  
  if (nrow(vi_mv_raw) == 0) {
    message("  No data for ", mv, " — skipping."); next
  }
  
  vi_mv_filtered <- filter_predictors(vi_mv_raw)
  
  pred_order_mv <- vi_mv_filtered %>%
    group_by(predictor) %>%
    summarise(med = median(importance_rel, na.rm = TRUE), .groups = "drop") %>%
    arrange(med) %>%
    pull(predictor)
  
  p2 <- make_vi_boxplot(
    df         = vi_mv_filtered,
    pred_order = pred_order_mv,
    subtitle   = paste0("BASE + ", MGMT_COL_LABELS[mv])
  )
  
  n_preds_mv <- length(pred_order_mv)
  save_fig(p2, name = paste0("plot2_", mv, "_RF"),  w = 30, h = 17)
  
}

# =============================================================================
# SECTION 6 — PI importance report (console summary)
# =============================================================================
# Prints mean relative importance of PI [0, 1] across gap sizes for:
#   (A) The managed run   — PI stored as PI_{gap_label} columns
#   (B) The phytomass_index ablation run — PI stored as PI_{gap_label}

message("\n", strrep("=", 60))
message("PI (Phytomass Index) — mean relative importance summary")
message(strrep("=", 60))

read_pi_from_vi <- function(vi_root, site, source_label) {
  map_dfr(GAP_SIZES, function(gs) {
    path <- file.path(vi_root, "variable_importance",
                      glue("rf_variable_importance_{gs}_all.rds"))
    if (!file.exists(path)) {
      message("  Not found: ", path); return(NULL)
    }
    raw <- readRDS(path)
    
    if (gs == GAP_SIZES[1]) {
      pi_like <- grep("^PI|phytomass", unique(raw$predictor),
                      value = TRUE, ignore.case = FALSE)
      message("  [", source_label, " | ", site, " | ", gs, "] PI-like: ",
              if (length(pi_like)) paste(head(pi_like, 10), collapse = ", ")
              else "(none found)")
    }
    
    pi_rows <- raw %>%
      filter(str_starts(predictor, "PI") |
               predictor == "PI"           |
               predictor == "phytomass_index")
    
    if (nrow(pi_rows) == 0) return(NULL)
    
    pi_rows %>%
      group_by(gap_size) %>%
      summarise(mean_imp = mean(importance_rel, na.rm = TRUE), .groups = "drop") %>%
      mutate(site = site, source = source_label, gap_size = gs)
  })
}

pi_managed <- map_dfr(SITES, function(s)
  read_pi_from_vi(file.path(RESULTS_ROOT, s, "managed", "RF"),
                  s, "managed"))

pi_effect <- map_dfr(SITES, function(s)
  read_pi_from_vi(file.path(RESULTS_ROOT, s, "management_effect",
                            "RF", "phytomass_index"),
                  s, "management_effect/phytomass_index"))

pi_report <- bind_rows(pi_managed, pi_effect)
message("PI report rows: ", nrow(pi_report))

if (nrow(pi_report) == 0) {
  message("No PI data found.")
} else {
  pi_out <- pi_report %>%
    mutate(
      mean_imp_fmt = sprintf("%.4f", mean_imp),   # [0, 1] — no * 100
      gap_size     = factor(gap_size, levels = GAP_SIZES)
    ) %>%
    select(source, site, gap_size, mean_imp_fmt) %>%
    arrange(source, site, gap_size)
  
  message(sprintf("%-45s %-5s %-8s %s",
                  "source", "site", "gap_size", "mean_imp_rel"))
  message(strrep("-", 70))
  for (i in seq_len(nrow(pi_out))) {
    r <- pi_out[i, ]
    message(sprintf("%-45s %-5s %-8s %s",
                    as.character(r$source),   as.character(r$site),
                    as.character(r$gap_size), r$mean_imp_fmt))
  }
}

message(strrep("=", 60))
message("\nAll done. Output: ", OUT_DIR)

# =============================================================================
# SECTION 7 — CSV importance summaries (median, Q25, Q75 per predictor)
# =============================================================================
# Produces 7 CSV files:
#   vi_summary_managed.csv          — managed predictor set
#   vi_summary_{MGMT_VAR}.csv       — one per feature-addition run
#
# Each CSV has columns:
#   predictor, site, median, Q25, Q75

OUT_CSV <- here::here("graphs", "VI", "csv")
dir.create(OUT_CSV, recursive = TRUE, showWarnings = FALSE)

compute_vi_summary <- function(df) {
  df %>%
    group_by(predictor, site) %>%
    summarise(
      median = median(importance_rel, na.rm = TRUE),
      Q25    = quantile(importance_rel, 0.25, na.rm = TRUE),
      Q75    = quantile(importance_rel, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(site, desc(median))
}

# --- Managed ---
message("\n[CSV] Writing managed summary ...")
compute_vi_summary(vi_managed_filtered) %>%
  write.csv(file.path(OUT_CSV, "vi_summary_managed.csv"), row.names = FALSE)
message("  Saved: vi_summary_managed.csv")

# --- Feature-addition runs ---
message("\n[CSV] Writing feature-addition summaries ...")

for (mv in MGMT_VARS) {
  
  vi_mv_raw <- map_dfr(SITES, function(s) {
    vi_root <- file.path(RESULTS_ROOT, s, "management_effect", "RF", mv)
    raw     <- load_rf_vi_raw(vi_root)
    if (is.null(raw) || nrow(raw) == 0) return(NULL)
    tidy_vi_raw(raw) %>% mutate(site = s)
  })
  
  if (nrow(vi_mv_raw) == 0) {
    message("  No data for ", mv, " — skipping."); next
  }
  
  fname <- paste0("vi_summary_", mv, ".csv")
  filter_predictors(vi_mv_raw) %>%
    compute_vi_summary() %>%
    write.csv(file.path(OUT_CSV, fname), row.names = FALSE)
  message("  Saved: ", fname)
}

message("\nCSV summaries written to: ", OUT_CSV)

# ====================== end plot_vi_summary.R ================================

# ====================== end plot_vi_summary.R ================================