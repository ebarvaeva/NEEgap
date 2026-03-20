# NEEgap

A complete R pipeline for gap-filling eddy-covariance NEE (Net Ecosystem Exchange) at managed grassland sites, comparing machine-learning models against process-based benchmarks and quantifying the contribution of management variables (grazing, fertilisation, canopy structure) to gap-filling accuracy.

Developed for sites **JC1** and **JC2** (Johnstown Castle, Teagasc, Ireland), 2023–2024.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Repository Structure](#2-repository-structure)
3. [Data Layout](#3-data-layout)
4. [Prerequisites](#4-prerequisites)
5. [Quick Start](#5-quick-start)
6. [Pipeline in Detail](#6-pipeline-in-detail)
   - [Stage 0 — Quality Control](#stage-0--quality-control)
   - [Stage 1 — Data Preparation](#stage-1--data-preparation)
   - [Stage 2 — Gap-Filling Models](#stage-2--gap-filling-models)
   - [Stage 3 — Management Ablation Study](#stage-3--management-ablation-study)
   - [Stage 4 — Metrics and Plots](#stage-4--metrics-and-plots)
7. [Models](#7-models)
8. [Cross-Validation Design](#8-cross-validation-design)
9. [Predictors](#9-predictors)
10. [Output Structure](#10-output-structure)
11. [Adding a New Site](#11-adding-a-new-site)

---

## 1. Overview

```
Raw EddyPro / Biomet / FluxNet CSVs
             │
             ▼
   [0]  Quality Control          quality_control_nee.R
             │  merged_qc.csv
             ▼
   [1]  Data Preparation         01–07_*.R
             │  data/JC1.rds,  JC2.rds
             ▼
   [2]  Gap-Filling Models       run_model.R  /  run_MDS_miniRECgap.R
     ┌───────┬───────┬───────┬────────────┬──────┐
     RF     MLP    XGB   miniRECgap    MDS
             │
             ▼
   [3]  Management Ablation      run_management_effect_*.R
        (BASE + each management variable, in isolation)
             │
             ▼
   [4]  Metrics & Plots          metrics/*.R
        (MAE, RMSE, R²  ·  variable importance  ·  management effect)
```

**What the pipeline does:**

- Fills gaps in half-hourly NEE measurements using five independent models.
- Derives Reco and GPP post-hoc from the gap-filled NEE through flux partitioning, guaranteeing the carbon-balance identity NEE = Reco − GPP at every gap row.
- Evaluates each model at four gap sizes (S ≈ 3 d, M ≈ 7 d, L ≈ 14 d, VL ≈ 30 d) via a leave-one-gap-out cross-validation.
- Isolates the contribution of each management variable to accuracy using a feature-addition ablation study.

---

## 2. Repository Structure

```
NEEgap/
│
├── data/                               # Raw and prepared data (not in repo)
│   ├── raw_data/                       # EddyPro, biomet, FluxNet, meta CSVs
│   ├── data_qc/                        # QC output (merged_qc.csv per site-year)
│   ├── data_prepared/                  # Final RDS files (JC1.rds, JC2.rds)
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
│   ├── 01_run_data_preparation.R       # [1] Master pipeline runner
│   ├── 02_timestamp_correction.R       # Full-year 30-min grid, NA rows for gaps
│   ├── 03_predictors_meteo.R           # PPFD, Rg, Temp, RH, VPD, rain, time encodings
│   ├── 04_management_days_since.R      # Grazing_days_since, Fertiliser_days_since
│   ├── 05_nitrogen.R                   # N step function from fertiliser events
│   ├── 06_grass_canopy.R               # grass_height, grass_biomass (interpolated)
│   ├── 07_bind_years.R                 # Bind 2023 + 2024, validate, save final RDS
│   ├── ManagementEvents.R              # Hard-coded JC1/JC2 management event tables
│   ├── generate_management_csvs.R      # One-time: export event tables to CSV files
│   └── Regrowth_period.R               # Regrowth period definition utilities
│
├── models/
│   ├── RF_CV.R                         # Random Forest cross-validation
│   ├── MLP_CV.R                        # Multi-Layer Perceptron cross-validation
│   ├── XGBoost_CV.R                    # XGBoost cross-validation
│   ├── miniRECgap_CV.R                 # Process-based Reco/GPP model
│   ├── MDS_CV.R                        # Marginal Distribution Sampling (REddyProc)
│   ├── RF_Grazing_CV.R                 # RF: BASE + Grazing_days_since, PI disabled
│   ├── RF_PI_CV.R                      # RF: BASE + Phytomass Index only
│   ├── MLP_Grazing_CV.R                # MLP ablation variants (same logic)
│   ├── MLP_PI_CV.R
│   ├── XGBoost_Grazing_CV.R
│   └── XGBoost_PI_CV.R
│
├── scripts/
│   ├── run_model.R                     # [2] Entry point: RF | MLP | XGBoost
│   ├── run_MDS_miniRECgap.R            # [2] Entry point: MDS | miniRECgap
│   ├── run_management_effect_RF.R      # [3] RF ablation: BASE + each management var
│   ├── run_management_effect_MLP.R     # [3] MLP ablation
│   ├── run_management_effect_XGBoost.R
│   ├── run_management_effect_GrazingDaysSince.R  # Grazing as direct feature, no PI
│   └── run_management_effect_PI.R      # PI only, no Grazing_days_since in predictors
│
├── metrics/
│   ├── compute_metrics_and_plots.R     # [4a] MAE/RMSE/R², overall and site comparison
│   ├── graph_colour_individual_gaps.R  # [4b] Season- and recovery-encoded plots
│   ├── management_effect_graphs.R      # [4c] Delta metrics, heatmaps, bump charts
│   └── variable_importance_graphs.R    # [4d] RF variable importance figures
│
├── graphs/                             # All plot outputs (created automatically)
└── results/                            # All model prediction outputs (created automatically)
```

---

## 3. Data Layout

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

## 4. Prerequisites

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

## 5. Quick Start

Open R with the project root as the working directory (e.g. open `GapFillNEE.Rproj` in RStudio).

```r
# One-time: generate management event CSV files
source("data_preparation/generate_management_csvs.R")

# 0 — Run QC (once per site × year, with paths set inside the script)
source("data_preparation/quality_control_nee/quality_control_nee.R")

# 1 — Prepare both sites, both years
source("data_preparation/01_run_data_preparation.R")

# 2a — Gap-fill with Random Forest, JC1, managed predictor set
#      Edit SITE_NAME, MANAGEMENT_CONDITION, MODEL_CHOICE in run_model.R first
source("scripts/run_model.R")

# 2b — Gap-fill with miniRECgap or MDS
#      Edit SITE_NAME and MODEL_CHOICE in run_MDS_miniRECgap.R first
source("scripts/run_MDS_miniRECgap.R")

# 3 — Management variable ablation study
source("scripts/run_management_effect_RF.R")
source("scripts/run_management_effect_MLP.R")
source("scripts/run_management_effect_XGBoost.R")
source("scripts/run_management_effect_GrazingDaysSince.R")  # Grazing only, no PI
source("scripts/run_management_effect_PI.R")                # PI only

# 4 — Metrics and plots (run in this order)
source("metrics/compute_metrics_and_plots.R")    # must run first — produces CSVs
source("metrics/graph_colour_individual_gaps.R")
source("metrics/management_effect_graphs.R")
source("metrics/variable_importance_graphs.R")
```

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

---

### Stage 2 — Gap-Filling Models

**Entry point:** `scripts/run_model.R` (RF, MLP, XGBoost) or `scripts/run_MDS_miniRECgap.R` (MDS, miniRECgap)

Edit the `USER SETTINGS` block at the top of the chosen script, then source it:

```r
# run_model.R — edit these three lines:
SITE_NAME            <- "JC1"       # "JC1" or "JC2"
MANAGEMENT_CONDITION <- "managed"   # "managed" or "unmanaged"
MODEL_CHOICE         <- "RF"        # "RF" | "MLP" | "XGBoost"
```

The script loads `data/{SITE_NAME}.rds`, intersects `FEATURE_SET` with the available columns, creates the output directory, and sources the appropriate model CV script in an isolated environment, passing `df`, `predictors`, `RESULTS_DIR`, and `rds_name`.

**Output:** `results/{SITE}/{MANAGEMENT}/{MODEL}/df_cv_all_predictions.rds`

---

### Stage 3 — Management Ablation Study

**Scripts:** `scripts/run_management_effect_{RF|MLP|XGBoost}.R`

Trains each model with the **BASE predictor set** and then, separately, with **BASE + one management variable at a time**, across both sites. This produces an isolated estimate of each variable's marginal contribution.

| Script | Predictor set | PI activated? |
|---|---|---|
| `run_management_effect_{RF|MLP|XGBoost}.R` | BASE + {Grazing\_days\_since, Fertiliser\_days\_since, N, grass\_height, grass\_biomass} each in turn | Auto (yes when Grazing is added) |
| `run_management_effect_GrazingDaysSince.R` | BASE + Grazing\_days\_since | **No** — Grazing as direct numeric feature only |
| `run_management_effect_PI.R` | BASE only | **Yes** — PI computed per gap, Grazing not in predictors |

**Delta convention:** `delta = metric(BASE + mgmt) − metric(BASE)` — negative delta for MAE/RMSE means improvement; positive for R².

**Output:** `results/{SITE}/management_effect/{MODEL}/{MGMT_VAR}/df_cv_all_predictions.rds`

---

### Stage 4 — Metrics and Plots

Run the scripts in this order — each depends on outputs from the previous:

```
compute_metrics_and_plots.R   →  produces metrics_csv/ used by the other three
graph_colour_individual_gaps.R
management_effect_graphs.R
variable_importance_graphs.R
```

| Script | Plots produced |
|---|---|
| `compute_metrics_and_plots.R` | Line-profile overall metrics (x = gap size, colour = model); connected-dot site-comparison plots |
| `graph_colour_individual_gaps.R` | Same connected-dot plots with line colour = dominant season, dot shape = grazing recovery phase |
| `management_effect_graphs.R` | ΔMAE/RMSE/R² barplots; absolute metric connected-dot plots; % MAE-reduction heatmaps; recovery-phase interaction; seasonal breakdown; ranking stability bump charts |
| `variable_importance_graphs.R` | RF impurity importance: ranked barplots; heatmaps; rank-position distributions; management share of total importance; per-gap timelines; cross-site summary |

All plots are saved as both **PNG** (no title, for figures) and **PDF** (with title, for review) under `graphs/`.

---

## 7. Models

| Model | Type | Package | Key details |
|---|---|---|---|
| **RF** | Random Forest | `ranger` | Impurity importance saved per gap label; no feature scaling required |
| **MLP** | Neural Network | `keras3` / TensorFlow | 128–64–32 ReLU layers; Adam optimiser; early stopping patience = 10 on validation loss; min-max scaling of features and target; stratified 10 % validation split by season × month × day/night |
| **XGBoost** | Gradient Boosting | `xgboost` | Early stopping on validation RMSE; boosting evaluation curves saved per gap label |
| **miniRECgap** | Process-based | Base R | Lloyd-Taylor Reco + Thornley non-rectangular hyperbola GPP; parameters fitted per regrowth period from outside-gap observations; produces Reco and GPP directly (not post-hoc) |
| **MDS** | Look-up table | `REddyProc` | Marginal Distribution Sampling (Wutzler et al., 2018); window ±5/10/20 days; requires a complete annual timeline; excluded from L and VL gap sizes |

### Reco and GPP derivation

For RF, MLP, and XGBoost, Reco and GPP are derived **post-hoc** from the gap-filled NEE through flux partitioning:

- **Night-time** (PPFD < 10 µmol m⁻² s⁻¹): Reco fitted by OLS to the Lloyd-Taylor temperature response; base respiration R10 estimated from nighttime training rows.
- **Daytime**: GPP fitted by BFGS to the Thornley non-rectangular hyperbola using daytime training rows.
- Parameters are estimated **per regrowth period** (delineated by grazing events) from observations outside each artificial gap, preserving cross-validation integrity and ensuring NEE = Reco − GPP holds at every gap row.

---

## 8. Cross-Validation Design

Each year-long dataset is partitioned into **contiguous artificial gaps** of four target durations:

| Label | Duration | Tolerance |
|---|---|---|
| S (Short) | ≈ 3 days (144 half-hours) | ± 10 rows |
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
| `Rg` | Global shortwave radiation | W m⁻² |
| `VPD` | Vapour pressure deficit | kPa |
| `RH` | Relative humidity | % |
| `Temp` | Air temperature | °C |
| `rain` | Precipitation per half-hour | mm |
| `rain_rolling_24` | 24 h trailing mean precipitation | mm |
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

Including `Grazing_days_since` in the predictor list also activates the **Phytomass Index** — a per-gap derived predictor computed from the rolling 21-day balance of night-time and daytime NEE observations outside each artificial gap. PI ∈ [0, 1] reflects canopy regrowth state without leaking any information from inside the gap being filled.

To use `Grazing_days_since` as a raw predictor *without* PI, use `run_management_effect_GrazingDaysSince.R`.
To use PI *without* `Grazing_days_since` in the predictor set, use `run_management_effect_PI.R`.

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
    └── management_effect/{MODEL}/{MGMT_VAR}/

graphs/
├── metrics_csv/                        # Per-gap metrics CSVs (input to other scripts)
│   └── gap_metrics_{TARGET}_{SIZE}.csv
├── overall_metrics/{TARGET}/           # Line-profile plots by gap size
├── site_comparison/{TARGET}/{SIZE}/    # Connected-dot comparison plots
├── site_comparison_coloured/{TARGET}/{SIZE}/  # Season/recovery-encoded version
├── management_effect/
│   ├── delta_barplots/
│   ├── connected_dots/
│   ├── heatmaps/
│   ├── recovery_interaction/
│   ├── seasonal_breakdown/
│   └── ranking_stability/
└── variable_importance/
    ├── 00_managed/,  00_unmanaged/
    ├── 01_ranked_barplots/
    ├── 02_heatmaps/
    ├── 03_rank_distributions/
    ├── 04_mgmt_share/
    ├── 05_timelines/
    └── 06_crosssite_summary/
```

---

## 11. Adding a New Site

1. **Place raw data** in `data/raw_data/{NEWSITE}_{YEAR}_*.csv` following the layout in [Section 3](#3-data-layout).

2. **QC:** in `quality_control_nee.R`, set `site_key` to the site's row name in `nasco_site_list_wkt.csv` and update the four `*_path` variables. Run for each year.

3. **Management events:** add the site's event tables to `ManagementEvents.R`, then re-run `generate_management_csvs.R`. Set `HAS_MANAGEMENT <- FALSE` in `01_run_data_preparation.R` to skip scripts 04–06 for an unmanaged site.

4. **Data preparation:** add the new site to `SITES` in `01_run_data_preparation.R`:
   ```r
   SITES <- c("JC1", "JC2", "NEWSITE")
   ```

5. **Gap-filling:** set `SITE_NAME <- "NEWSITE"` in `run_model.R` and run for each model and management condition.

6. **Ablation study:** add `"NEWSITE"` to the `SITES` vector in each `run_management_effect_*.R` script.

7. **Metrics:** add the new site to `SITES_ALL` in `compute_metrics_and_plots.R` and the other three metrics scripts.

---

## References

- Kljun, N., Calanca, P., Rotach, M.W., Schmid, H.P. (2015). A simple two-dimensional parameterisation for Flux Footprint Prediction (FFP). *Geoscientific Model Development*, 8, 3695–3713. https://doi.org/10.5194/gmd-8-3695-2015
- McCree, K.J. (1972). The action spectrum, absorptance and quantum yield of photosynthesis in crop plants. *Agricultural Meteorology*, 9, 191–216.
- Wutzler, T., Lucas-Moffat, A., Migliavacca, M., Knauer, J., Sickel, K., Šigut, L., Menzer, O., Reichstein, M. (2018). Basic and extensible post-processing of eddy covariance flux data with EddyPro in combination with REddyProc. *Biogeosciences*, 15, 5015–5030. https://doi.org/10.5194/bg-15-5015-2018
