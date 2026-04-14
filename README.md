# NEEgap

A complete R pipeline for gap-filling eddy-covariance NEE (Net Ecosystem Exchange) at managed grassland sites, comparing machine-learning models against process-based benchmarks and quantifying the contribution of management variables (grazing, fertilisation, canopy structure) to gap-filling accuracy.

Developed for sites **JC1** and **JC2** (Johnstown Castle, Teagasc, Ireland), 2023–2024.

---

## Table of Contents

1. [Quick Start](#1-quick-start)
2. [Overview](#2-overview)
3. [Repository Structure](#3-repository-structure)
4. [Data Layout](#4-data-layout)
5. [Prerequisites](#5-prerequisites)
6. [Pipeline in Detail](#6-pipeline-in-detail)
   - [Stage 0 — Quality Control](#stage-0--quality-control)
   - [Stage 1 — Data Preparation](#stage-1--data-preparation)
   - [Stage 2 — Gap-Filling Models](#stage-2--gap-filling-models)
   - [Stage 3 — Management Effect Study](#stage-3--management-effect-study)
   - [Stage 4 — Metrics and Plots](#stage-4--metrics-and-plots)
7. [Models](#7-models)
8. [Cross-Validation Design](#8-cross-validation-design)
9. [Predictors](#9-predictors)
10. [Output Structure](#10-output-structure)
11. [Adding a New Site](#11-adding-a-new-site)

---

## 1. Quick Start

Open R with the project root as the working directory (e.g. open `GapFillNEE.Rproj` in RStudio).

```r
# One-time: generate management event CSV files
source("data_preparation/generate_management_csvs.R")

# 0 — Run QC (once per site × year, with paths set inside the script)
source("data_preparation/quality_control_nee/quality_control_nee.R")

# 1a — Prepare both sites, both years → JC1.rds, JC2.rds
source("data_preparation/01_run_data_preparation.R")

# 1b — Construct artificial gap flags and masked-NEE columns → JC1_cv.rds, JC2_cv.rds
source("data_preparation/08__artificial_gaps.R")

# 1c — Compute Phytomass Index per gap label → appended to JC1_cv.rds, JC2_cv.rds
source("data_preparation/09__phytomass_index.R")

# 2a — Gap-fill with Random Forest, JC1, managed predictor set
#      Edit SITE_NAME, MANAGEMENT_CONDITION, MODEL_CHOICE in run_model.R first
source("scripts/run_model.R")

# 2b — Gap-fill with miniRECgap or MDS
#      Edit SITE_NAME and MODEL_CHOICE in run_MDS_miniRECgap.R first
source("scripts/run_MDS_miniRECgap.R")

# 3 — Management variable feature addition study
#     Edit SITE_NAME, MODEL_CHOICE, and MGMT_EVENT in run_management_effect.R,
#     then source once per combination.
#     MGMT_EVENT options: "Grazing_days_since" | "Fertiliser_days_since" |
#                         "N" | "grass_height" | "grass_biomass" | "PI"
source("scripts/run_management_effect.R")

# 4 — Metrics and plots
# Run compute_metrics_and_plots.R first — it produces the CSVs used by
# graph_colour_individual_gaps.R and management_effect_graphs.R.
# The remaining scripts are independent and can be run in any order.
source("metrics/compute_metrics_and_plots.R")       # must run first — produces CSVs
source("metrics/graph_colour_individual_gaps.R")
source("metrics/management_effect_graphs.R")
source("metrics/VI_graphs_compact.R")
source("metrics/PI_threshold_plot.R")
source("metrics/plot_artificial_gaps.R")
source("metrics/plot_real_gaps_cleaveland.R")
source("metrics/timeseries_nee.R")
```

---

## 2. Overview

```
Raw EddyPro / Biomet / FluxNet CSVs
             │
             ▼
   [0]  Quality Control          quality_control_nee.R
             │  merged_qc.csv
             ▼
   [1]  Data Preparation         01–09_*.R
             │  data/JC1_cv.rds,  JC2_cv.rds
             ▼
   [2]  Gap-Filling Models       run_model.R  /  run_MDS_miniRECgap.R
     ┌───────┬───────┬───────┬────────────┬──────┐
     RF     MLP    XGB   miniRECgap    MDS
             │
             ▼
   [3]  Management Effect (Feature Addition)      run_management_effect.R
        (BASE + each management variable)
             │
             ▼
   [4]  Metrics & Plots          metrics/*.R
        (MAE, RMSE, R²  ·  variable importance  ·  management effect)
```

**What the pipeline does:**

- Fills gaps in half-hourly NEE measurements using five independent models.
- Evaluates each model at four gap sizes (S ≈ 1 d, M ≈ 7 d, L ≈ 14 d, VL ≈ 30 d) via a leave-one-gap-out cross-validation.
- Isolates the contribution of each management variable to accuracy using a feature-addition ablation study.

---

## 3. Repository Structure

```
NEEgap/
│
├── data/                               # Raw and prepared data (not in repo)
│   ├── raw_data/                       # EddyPro, biomet, FluxNet, meta CSVs
│   ├── data_qc/                        # QC output (merged_qc.csv per site-year)
│   ├── data_prepared/                  # Prepared RDS files (JC1.rds, JC2.rds, JC1_cv.rds, JC2_cv.rds)
│   ├── management_data/                # Grazing/fertiliser event CSVs
│   ├── canopy_data/                    # Grass height and biomass measurements
│   ├── meteireann_data/                # Met Éireann hourly weather and solar
│   └── site_polygon/                   # Site boundary coordinates
│
├── data_preparation/
│   ├── quality_control_nee/
│   │   ├── quality_control_nee.R       # [0] QC pipeline (base filters + footprint)
│   │   └── grid_qc/
│   │       ├── pblh_calc.R             # Planetary boundary layer height (Kljun 2015)
│   │       ├── calc_footprint_FFP_mod.R  # 2-D FFP footprint model (Kljun 2015)
│   │       └── utilities.R             # Spatial grid helpers
│   ├── 01_run_data_preparation.R       # [1] Master pipeline runner (scripts 02–07)
│   ├── 02_timestamp_correction.R       # Full-year 30-min grid, NA rows for gaps
│   ├── 03_predictors_meteo.R           # PPFD, Rg, Temp, RH, VPD, rain, time encodings
│   ├── 04_management_days_since.R      # Grazing_days_since, Fertiliser_days_since
│   ├── 05_nitrogen.R                   # N step function from fertiliser events
│   ├── 06_grass_canopy.R               # grass_height, grass_biomass (interpolated)
│   ├── 07_bind_years.R                 # Bind 2023 + 2024, validate, save final RDS
│   ├── 08__artificial_gaps.R           # [1b] Construct CV gap flags and masked-NEE columns → JCi_cv.rds
│   ├── 09__phytomass_index.R           # [1c] Compute PI_{gap_label} columns → appended to JCi_cv.rds
│   ├── ManagementEvents.R              # Hard-coded JC1/JC2 management event tables
│   ├── generate_management_csvs.R      # One-time: export event tables to CSV files
│   └── Regrowth_period.R               # Regrowth period definition utilities
│
├── models/
│   ├── RF_CV.R                         # Random Forest cross-validation
│   ├── MLP_CV.R                        # Multi-Layer Perceptron cross-validation
│   ├── XGBoost_CV.R                    # XGBoost cross-validation
│   ├── miniRECgap_CV.R                 # Process-based Reco/GPP model
│   └── MDS_CV.R                        # Marginal Distribution Sampling (REddyProc)
│
├── scripts/
│   ├── run_model.R                     # [2] Entry point: RF | MLP | XGBoost
│   ├── run_MDS_miniRECgap.R            # [2] Entry point: MDS | miniRECgap
│   └── run_management_effect.R         # [3] Feature Addition: set MODEL_CHOICE + MGMT_EVENT
│
├── metrics/
│   ├── compute_metrics_and_plots.R     # [4a] MAE/RMSE/R² (pooled & split), CSV export
│   ├── graph_colour_individual_gaps.R  # [4b] Season- & recovery-encoded connected-dot plots
│   ├── management_effect_graphs.R      # [4c] Management-effect CSVs + % MAE-reduction heatmaps
│   ├── VI_graphs_compact.R             # [4d] RF variable importance boxplot summaries
│   ├── PI_threshold_plot.R             # [4e] PI sensitivity: PPFD > 400 vs > 700 scatter
│   ├── plot_artificial_gaps.R          # [4f] Artificial gap positions on NEE time series
│   ├── plot_real_gaps_cleaveland.R     # [4g] Real gap length distribution (Cleveland dot plot, log scale)
│   └── timeseries_nee.R                # [4h] NEE time series + management events + scatter
│
├── graphs/                             # All plot outputs (created automatically)
└── results/                            # All model prediction outputs (created automatically)
```

---

## 4. Data Layout

Data files are **not included** in this repository. Place them in the `data/` directory:

```
data/
├── raw_data/
│   ├── {SITE}_{YEAR}_eddypro.csv      # EddyPro full output
│   ├── {SITE}_{YEAR}_biomet.csv       # Biomet half-hourly
│   ├── {SITE}_{YEAR}_fluxnet.csv      # FluxNet-format output
│   └── {SITE}_{YEAR}_meta.csv         # Site metadata (sonic height, roughness length, lat)
│
├── meteireann_data/
│   ├── met_hourly.csv                 # Hourly Temp, RH, rain  (date DD/MM/YYYY HH:MM)
│   └── solar_hourly.xlsx              # Hourly global radiation  (date, glorad W m⁻²)
│
├── management_data/                   # Generated by generate_management_csvs.R
│   ├── grazing_events_{SITE}.csv      # date, start_time, end_time, event_id
│   ├── fertiliser_events_{SITE}.csv
│   ├── last_events_before_year_{SITE}.csv  # year_for, event_type, last_end
│   └── nitrogen_amounts_{SITE}.csv    # fert_date, N_kg_ha, is_slurry
│
├── canopy_data/
│   └── grass_drybiomass_JC_2023_2024.rds  # date, height_JC1/JC2, drymass_JC1/JC2
│
└── site_polygon/
    └── {SITE}_polygon.csv             # lon, lat (closed ring, for footprint QC)
```

**Generate management event CSVs once** from the hard-coded event tables in `ManagementEvents.R`:

```r
source("data_preparation/generate_management_csvs.R")
```

After this, CSVs can be edited manually for new years or sites.

---

## 5. Prerequisites

**R ≥ 4.2** with the following packages:

| Scope | Packages |
|---|---|
| Core (all scripts) | `here`, `dplyr`, `tibble`, `tidyr`, `purrr`, `lubridate`, `readr`, `stringr`, `ggplot2`, `zoo` |
| Random Forest | `ranger` |
| MLP | `reticulate`, `tensorflow`, `keras3` |
| XGBoost | `xgboost` |
| MDS | `REddyProc` |
| QC footprint | `sf` |
| Data reading | `readxl` |
| Plotting | `patchwork` |

**Python / TensorFlow** (MLP only): create a virtual environment named `r-tensorflow` with TensorFlow ≥ 2.x, then configure `reticulate` to use it. MLP can be skipped without affecting any other step.

---

## 6. Pipeline in Detail

### Stage 0 — Quality Control

**Script:** `data_preparation/quality_control_nee/quality_control_nee.R`

Reads four EddyPro-generated input files (eddypro, biomet, meta, fluxnet), merges them onto a single half-hourly timeline, and applies a multi-stage QC procedure.

**Base filters** (each produces a separate diagnostic column for inspection):

| Filter | Criterion |
|---|---|
| CO₂ range | `co2_flux` outside [−40, +30] µmol m⁻² s⁻¹ |
| Signal strength | LI-7200 CO₂ signal strength < 70 % |
| Random error | `rand_err_co2_flux` > 100 µmol m⁻² s⁻¹ |
| Flow rate | Sample tube flow outside [12, 18] L min⁻¹ |

**Stationarity / turbulence flag:** EddyPro's `qc_co2_flux` column (0–9 scale, Foken et al. 2004) is applied at all nine threshold levels, producing columns `co2_flux_base_filters_{1..9}`.

**Grid footprint QC** (Kljun et al., 2015): PBLH is computed per half-hour (`pblh_calc.R`), then the 2-D flux footprint probability density is generated (`calc_footprint_FFP_mod.R`). The fraction of footprint mass inside the site boundary polygon is computed with `sf::st_within()`, and three threshold flags (70 / 80 / 90 %) are applied, producing columns `co2_flux_base_filters_{i}_{p}_grid`.

**The data preparation pipeline reads `co2_flux_base_filters_6_70_grid`** as the standard analysis column (stationarity flag ≤ 6, footprint ≥ 70 %).

**Outputs:** `data/data_qc/results_qc_{SITE}_{YEAR}/merged_qc.csv` and `merged_qc.rds`

---

### Stage 1 — Data Preparation

**Stage 1a — Predictor construction**

**Master script:** `data_preparation/01_run_data_preparation.R`

Sources sub-scripts 02–07 in sequence for each site × year. Set `HAS_MANAGEMENT = FALSE` in the configuration block for unmanaged sites to skip scripts 04–06.

| Script | What it builds |
|---|---|
| `02` | Complete 30-min calendar-year grid; NA rows for measurement gaps |
| `03` | PPFD (3-priority fill chain), Rg, Temp, RH, VPD, rain, rain\_rolling\_24, night flag, temporal sine/cosine encodings, season dummies |
| `04` | `Grazing_days_since`, `Fertiliser_days_since`, `Fertiliser` dummy — continuous, 0 during events, back-filled across the year boundary |
| `05` | `N` — step function of kg N ha⁻¹ at each mineral fertiliser application |
| `06` | `grass_height` (cm), `grass_biomass` (kg DM ha⁻¹) — linearly interpolated from sparse field measurements, with a post-grazing anchor of 4 cm |
| `07` | Concatenates 2023 + 2024, validates zero NA in all non-NEE predictors, saves `data/data_prepared/{SITE}.rds` |

**PPFD gap-filling priority:**

1. On-site biomet (`PPFD_1_1_1`)
2. Twin-site biomet (the other grassland site, ≈ 2 km away)
3. Met Éireann global radiation: PPFD = Rg × 0.45 × 4.57  (McCree, 1972)

**Stage 1b — Artificial gap construction**

**Script:** `data_preparation/08__artificial_gaps.R`

Reads `{SITE}.rds`, filters to rows with finite `NEE_orig`, and constructs four sets of contiguous leave-one-gap-out (LOGO) blocks. For each gap label (e.g. `M3`) two columns are added:

- `M3` — logical flag; `TRUE` where this gap falls
- `NEE_M3` — `NEE_orig` with gap rows set to `NA` (training target column)

**Output:** `data/data_prepared/{SITE}_cv.rds` (one file per site, all four gap-size sets combined)

**Stage 1c — Phytomass Index**

**Script:** `data_preparation/09__phytomass_index.R`

Reads `{SITE}_cv.rds` and appends one `PI_{gap_label}` column per gap label across all four gap sizes (VL / L / M / S). PI is computed exclusively from observations *outside* each gap to prevent information leakage, then linearly interpolated across gap rows.

**Output:** `data/data_prepared/{SITE}_cv.rds` — same file, overwritten with PI columns appended

---

### Stage 2 — Gap-Filling Models

**Entry point:** `scripts/run_model.R` (RF, MLP, XGBoost) or `scripts/run_MDS_miniRECgap.R` (MDS, miniRECgap)

Edit the `USER SETTINGS` block at the top of the chosen script, then source it:

```r
# run_model.R — edit these lines:
SITE_NAME            <- "JC1"       # "JC1" or "JC2"
MANAGEMENT_CONDITION <- "managed"   # "managed" or "unmanaged"
MODEL_CHOICE         <- "RF"        # "RF" | "MLP" | "XGBoost"
PI_ENABLED           <- TRUE        # TRUE to include PI_{gap_label} as a predictor
```

The script loads `data/data_prepared/{SITE_NAME}_cv.rds` (produced by scripts 08 and 09), which already contains the gap flag columns, masked-NEE columns, and PI columns. It intersects `FEATURE_SET` with the available columns, creates the output directory, and sources the appropriate model CV script in an isolated environment, passing `df`, `predictors`, `RESULTS_DIR`, `rds_name`, and `PI_ENABLED`.

**Output:** `results/{SITE}/{MANAGEMENT}/{MODEL}/df_cv_all_predictions.rds`

---

### Stage 3 — Management Effect Study

**Script:** `scripts/run_management_effect.R`

Trains the chosen model with the **BASE predictor set** and then with **BASE + one management variable**, producing an isolated estimate of each variable's marginal contribution. Run once per management event / model / site combination by setting the three variables below:

```r
# run_management_effect.R — edit these lines:
SITE_NAME    <- "JC1"                 # "JC1" or "JC2"
MODEL_CHOICE <- "RF"                  # "RF" | "MLP" | "XGBoost"
MGMT_EVENT   <- "Grazing_days_since"  # see table below
```

| `MGMT_EVENT` value | Predictor added to BASE | PI activated? |
|---|---|---|
| `"Grazing_days_since"` | `Grazing_days_since` | No |
| `"Fertiliser_days_since"` | `Fertiliser_days_since` | No |
| `"N"` | `N` | No |
| `"grass_height"` | `grass_height` | No |
| `"grass_biomass"` | `grass_biomass` | No |
| `"PI"` | *(none — BASE only)* | Yes |

**Delta convention:** `delta = metric(BASE + mgmt) − metric(BASE)` — negative delta for MAE/RMSE means improvement; positive for R².

**Output:** `results/{SITE}/management_effect/{MODEL}/{MGMT_EVENT}/df_cv_all_predictions.rds`

---

### Stage 4 — Metrics and Plots

Run `compute_metrics_and_plots.R` first — it writes the CSVs that `graph_colour_individual_gaps.R` and `management_effect_graphs.R` read. The remaining five scripts are independent.

```
compute_metrics_and_plots.R   →  produces metrics_csv/ used by:
    graph_colour_individual_gaps.R
    management_effect_graphs.R

VI_graphs_compact.R           }
PI_threshold_plot.R           }  independent — run in any order
plot_artificial_gaps.R        }
plot_real_gaps_cleaveland.R   }
timeseries_nee.R              }
```

| Script | Plots produced |
|---|---|
| `compute_metrics_and_plots.R` | Line-profile overall metrics (x = gap size, colour = model, rows = ≤30 d / >30 d grazing window); same plots pooled over all rows (no split); writes `gap_metrics_NEE_*.csv`, `gap_metrics_whole_NEE_*.csv`, and `overall_metrics_pooled_NEE.csv` |
| `graph_colour_individual_gaps.R` | Connected-dot per-gap plots with line colour = dominant season and dot shape = grazing-recovery fraction; two families: (A) split by ≤30 d / >30 d, (B) whole-gap unsplit |
| `management_effect_graphs.R` | Management-effect pooled CSVs per variable; % MAE-reduction heatmaps (variable × model, faceted by site), one PNG per gap size |
| `VI_graphs_compact.R` | RF variable importance boxplot summaries across sites and gap sizes; managed baseline and feature-addition runs; PI columns collapsed to unified `"PI"` label |
| `PI_threshold_plot.R` | PI sensitivity scatter: PPFD > 400 vs PPFD > 700 threshold comparison for JC1 and JC2, coloured by season |
| `plot_artificial_gaps.R` | Artificial gap positions overlaid on the NEE time series; per-site × gap-size PNGs and side-by-side JC1/JC2 comparison panels |
| `plot_real_gaps_cleaveland.R` | Cleveland dot plot of real NEE gap counts by length (0–30 d bins + >30 d) for JC1 vs JC2 on a log x-axis |
| `timeseries_nee.R` | (Fig 1) NEE time series with annotated management events per site; (Fig 2/3) JC1 vs JC2 NEE scatter plots coloured by days since fertiliser, faceted by grazing recovery window |

All plots are saved as **PNG** under `graphs/`; diagnostic figures also export **PDF**.

---

## 7. Models

| Model | Type | Package | Key details |
|---|---|---|---|
| **RF** | Random Forest | `ranger` | Impurity importance saved per gap label; no feature scaling required |
| **MLP** | Neural Network | `keras3` / TensorFlow | 128–64–32 ReLU layers; Adam optimiser; early stopping patience = 10 on validation loss; min-max scaling of features and target; stratified 10 % validation split by season × month × day/night |
| **XGBoost** | Gradient Boosting | `xgboost` | Early stopping on validation RMSE; boosting evaluation curves saved per gap label |
| **miniRECgap** | Process-based | Base R | Lloyd-Taylor Reco + Thornley non-rectangular hyperbola GPP; parameters fitted per regrowth period from outside-gap observations; produces Reco and GPP directly (not post-hoc) |
| **MDS** | Look-up table | `REddyProc` | Marginal Distribution Sampling; window ±5/10/20 days; requires a complete annual timeline; excluded from L and VL gap sizes |

All model CV scripts (`RF_CV.R`, `MLP_CV.R`, `XGBoost_CV.R`, `miniRECgap_CV.R`, `MDS_CV.R`) expect `df` to be the output of scripts 08 and 09 (`{SITE}_cv.rds`). Gap construction, NA filtering, and PI computation are handled upstream and are not repeated inside the model scripts.

---

## 8. Cross-Validation Design

Each year-long dataset is partitioned into **contiguous artificial gaps** of four target durations (by script 08):

| Label | Duration | Tolerance |
|---|---|---|
| S (Short) | ≈ 1 day (48 half-hours) | ± 10 rows |
| M (Medium) | ≈ 7 days (336 half-hours) | ± 60 rows |
| L (Long) | ≈ 14 days (672 half-hours) | ± 120 rows |
| VL (Very Long) | ≈ 30 days (1 440 half-hours) | ± 250 rows |

Gaps are non-overlapping and constructed sequentially (VL → L → M → S) so that all four gap matrices coexist in the same data frame. For each gap label the model is trained on all observed rows outside that gap and predicts the gap rows — a **leave-one-gap-out (LOGO)** cross-validation.

**Output prediction column naming:** `{NEE|Reco|GPP}_{S|M|L|VL}_{rf|mlp|xgb|minirec|mds}_predicted`

Metrics are reported separately for two temporal windows per gap label:

| Window | Definition |
|---|---|
| `le30` | Gap rows where `Grazing_days_since` ≤ 30 days (early post-grazing recovery) |
| `gt30` | Gap rows where `Grazing_days_since` > 30 days (late recovery / ungrazed periods) |

---

## 9. Predictors

### BASE set (meteorological + temporal)

| Column | Description | Units |
|---|---|---|
| `PPFD` | Photosynthetic photon flux density | µmol m⁻² s⁻¹ |
| `Rg` | Global radiation | W m⁻² |
| `VPD` | Vapour pressure deficit | kPa |
| `RH` | Relative humidity | % |
| `Temp` | Air temperature | °C |
| `rain` | Precipitation per half-hour | mm |
| `rain_rolling_24` | 24 h rolling mean precipitation | mm |
| `hour_sin`, `hour_cos` | Hour of day — cyclical encoding | — |
| `doy_sin`, `doy_cos` | Day of year — cyclical encoding | — |
| `month_sin`, `month_cos` | Month of year — cyclical encoding | — |
| `Winter`, `Spring`, `Summer`, `Autumn` | Season dummies (Irish meteorological calendar) | 0 / 1 |
| `night` | Night flag: 1 if PPFD < 10 µmol m⁻² s⁻¹ | 0 / 1 |

### Management predictors (added for managed sites)

| Column | Description | Units |
|---|---|---|
| `Grazing_days_since` | Days since last grazing event ended | days |
| `Fertiliser_days_since` | Days since last fertiliser event ended | days |
| `N` | N applied at the most recent event (step function) | kg N ha⁻¹ |
| `grass_height` | Interpolated sward height | cm |
| `grass_biomass` | Interpolated above-ground dry matter | kg DM ha⁻¹ |

### Phytomass Index (PI)

**Phytomass Index** — a per-gap derived predictor computed by script 09 from the rolling 21-day balance of night-time and daytime NEE observations outside each artificial gap. PI ∈ [0, 1] reflects canopy regrowth state without leaking any information from inside the gap being filled. Each model script selects the matching `PI_{gap_label}` column automatically (e.g. `PI_S1` for gap `S1`).

---

## 10. Output Structure

```
results/
└── {SITE}/
    ├── managed/{MODEL}/
    │   ├── df_cv_all_predictions.rds   # 12 prediction columns: NEE/Reco/GPP × S/M/L/VL
    │   ├── run_info.txt
    │   ├── progress.log
    │   └── variable_importance/        # RF only
    │       ├── rf_variable_importance_{SIZE}_all.{rds,csv}
    │       └── rf_vi_{SIZE}_{LABEL}.rds
    ├── unmanaged/{MODEL}/
    ├── miniRECgap/
    ├── MDS/
    └── management_effect/{MODEL}/{MGMT_EVENT}/

graphs/
├── metrics_csv/                             # CSVs written by compute_metrics_and_plots.R
│   ├── gap_metrics_NEE_{SIZE}.csv           # Per-gap metrics (≤30 d / >30 d split)
│   ├── gap_metrics_whole_NEE_{SIZE}.csv     # Per-gap metrics (whole gap, no split)
│   ├── overall_metrics_pooled_NEE.csv       # Pooled metrics across all gaps per-size (all models)
│   └── management_effect/{SITE}/{VAR}/      # Per-variable pooled CSVs
│       └── overall_metrics_pooled_NEE.csv
├── overall_metrics/NEE/                     # Line-profile plots (split + no-split)
├── site_comparison_coloured/NEE/{SIZE}/     # Season/recovery-encoded split plots
├── site_comparison_whole/NEE/{SIZE}/        # Whole-gap unsplit connected-dot plots
├── management_effect/
│   └── heatmaps/NEE/                        # % MAE-reduction heatmaps per gap size
├── variable_importance/                     # RF VI boxplots (VI_graphs_compact.R)
├── PI/                                      # PI threshold sensitivity scatter
│   └── threshold_comparison_JC1_JC2.png
├── artificial_gaps/                         # Gap-position overlays on NEE series
│   ├── {SITE}/gap_positions_{SITE}_{SIZE}.{png,pdf}
│   └── comparison/gap_positions_comparison_{SIZE}.{png,pdf}
├── real_gaps/histogram/                     # Cleveland dot plot of real gap lengths
│   └── real_gap_cleveland.{png,pdf}
└── nee_sites/                               # NEE time series + scatter (timeseries_nee.R)
```

---

## 11. Adding a New Site

1. **Place raw data** in `data/raw_data/{NEWSITE}_{YEAR}_*.csv` following the layout in [Section 4](#4-data-layout).

2. **QC:** in `quality_control_nee.R`, set `site_key` to the site's row name in `nasco_site_list_wkt.csv` and update the four `*_path` variables. Run for each year.

3. **Management events:** add the site's event tables to `ManagementEvents.R`, then re-run `generate_management_csvs.R`. Set `HAS_MANAGEMENT <- FALSE` in `01_run_data_preparation.R` to skip scripts 04–06 for an unmanaged site.

4. **Data preparation:** add the new site to `SITES` in `01_run_data_preparation.R`:
   ```r
   SITES <- c("JC1", "JC2", "NEWSITE")
   ```
   Then add it to `SITES` in `08__artificial_gaps.R` and `09__phytomass_index.R` as well.

5. **Gap-filling:** set `SITE_NAME <- "NEWSITE"` in `run_model.R` and run for each model and management condition.

6. **Management Effect study:** set `SITE_NAME <- "NEWSITE"` in `run_management_effect.R` and run for each model and management event combination.

7. **Metrics:** add the new site to `SITES_ALL` in `compute_metrics_and_plots.R`, `graph_colour_individual_gaps.R`, `management_effect_graphs.R`, and `VI_graphs_compact.R`. Update `SITES` in `plot_artificial_gaps.R`, `plot_real_gaps_cleaveland.R`, `PI_threshold_plot.R`, and `timeseries_nee.R` as well.

---

## References
