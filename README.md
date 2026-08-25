# NEEgap

Code for the study **"Incorporating Management and Canopy Information into Machine-Learning Gap-Filling of Eddy-Covariance Net Ecosystem Exchange in Irish Grasslands"** (Enola Barvaeva, Andrew Parnell, Rachael Murphy, Morad Mirzaei, Ryan Burger, Katarina Domijan).

It gap-fills half-hourly **NEE** at two rotationally grazed dairy-grassland eddy-covariance towers at Johnstown Castle — **JC1** (2020, 2023, 2024) and **JC2** (2023, 2024) — and tests whether adding **management and canopy predictors** (days since grazing / fertilisation, sward height, dry-matter biomass, a flux-derived phytomass index) to the usual meteorological plus temporal drivers improves reconstruction.

Five models are compared under a **leave-one-gap-out cross-validation** across four gap lengths (S ≈ 1 d, M ≈ 7 d, L ≈ 14 d, VL ≈ 30 d): Random Forest (RF), XGBoost, multilayer perceptron (MLP), and two process-based benchmarks, Marginal Distribution Sampling (MDS) from REddyProc and miniRECgap.

---

## Repository structure

```text
NEEgap/
├── data_preparation/
│   ├── quality_control_nee/
│   │   ├── JC1_{2020,2023,2024}_qc.R          per site-year flux screening
│   │   ├── JC2_{2023,2024}_qc.R
│   │   └── grid_qc/
│   │       ├── calc_footprint_FFP_mod.R       2-D footprint model (Kljun 2015)
│   │       ├── pblh_calc.R                    boundary-layer height
│   │       ├── utilities.R                    spatial-grid helpers
│   │       └── nasco_site_list_wkt.csv        field boundary polygons — not included (see Data availability)
│   ├── 01_run_data_preparation.R             master runner for 02–07
│   ├── 02_timestamp_correction.R             full 30-min grid, gaps → NA rows
│   ├── 03_predictors_meteo.R                 PPFD, Rg, Temp, RH, VPD, rain, time encodings
│   ├── 04_management_days_since.R            Grazing_days_since, Fertiliser_days_since
│   ├── 05_nitrogen.R                         N step function (kg N ha⁻¹)
│   ├── 06_grass_canopy.R                     grass_height, grass_biomass (interpolated)
│   ├── 07_bind_years.R                       bind years, validate → {SITE}.rds
│   ├── 08_artificial_gaps.R                  CV gap flags + masked NEE → {SITE}_cv.rds
│   ├── 09_phytomass_index.R                  per-gap PI columns → {SITE}_cv.rds
│   ├── ManagementEvents.R                    hard-coded event tables
│   ├── generate_management_csvs.R            export event tables → data/management_data/
│   └── Regrowth_period.R                     standalone utility (legacy paths; not in pipeline)
│
├── models/                                    one CV gap-filler each; called by run scripts
│   ├── RF_CV.R   MLP_CV.R   XGBoost_CV.R
│   └── miniRECgap_CV.R   MDS_CV.R
│
├── run_models/
│   ├── run_ML_models/          JC{1,2}_{rf,mlp,xgb}_{mng,nomng}.R       
│   ├── run_benchmark_models/   JC{1,2}_{MDS,miniRECgap}.R              
│   └── run_management_effect/  run_management_effect_{RF,MLP,XGBoost,PI}.R
│
└── figures/                                   paper figures + summary tables
    ├── metrics_graphs.R                       MAE/RMSE/R² + writes graphs/metrics_csv/  
    ├── metrics_for_individual_gaps_graphs.R   per-gap plots, season/grazing encoded
    ├── management_effect_graphs.R             % MAE-reduction heatmaps
    ├── VI_graphs.R                            RF variable-importance boxplots
    ├── real_gaps_cleaveland_graphs.R          real gap-length distribution (Cleveland dot plot)
    ├── timeseries_nee_graphs.R                NEE time series + JC1-vs-JC2 scatter
    ├── wind_rose_graphs.R                     wind roses per site × year
    └── montly_temperature_rainfall_table.R    monthly temp/rainfall table 
```

`results/` (model predictions) and `graphs/` are created automatically by the run and figure scripts.

---

## Quick start

Run from the repository root so `here::here()` resolves. Install the packages (see **Requirements**) and place the `data/` tree (see **Data layout**) first; Stage 0 also needs `grid_qc/nasco_site_list_wkt.csv`.

```r
# 0 — Quality control    raw EddyPro exports → data/data_qc/results_qc_{SITE}_{YEAR}/
source("data_preparation/quality_control_nee/JC1_2020_qc.R")
source("data_preparation/quality_control_nee/JC1_2023_qc.R")
source("data_preparation/quality_control_nee/JC1_2024_qc.R")
source("data_preparation/quality_control_nee/JC2_2023_qc.R")
source("data_preparation/quality_control_nee/JC2_2024_qc.R")

# 1 — Data preparation
source("data_preparation/01_run_data_preparation.R")   # runs 02–07 per site-year → {SITE}.rds
source("data_preparation/08_artificial_gaps.R")         # → {SITE}_cv.rds  (gap flags + masked NEE)
source("data_preparation/09_phytomass_index.R")         # → appends PI_{gap} columns to {SITE}_cv.rds

# 2 — Models             each writes results/{SITE}/.../df_cv_all_predictions.rds
#   2a — ML gap-fillers  run all 12 scripts (rf/mlp/xgb × mng/nomng × JC1/JC2)
source("run_models/run_ML_models/JC1_rf_mng.R")         # … and the other 11
#   2b — Benchmarks
source("run_models/run_benchmark_models/JC1_MDS.R")     # + JC1_miniRECgap.R, JC2_MDS.R, JC2_miniRECgap.R
#   2c — Management effect (BASE + one variable at a time; each loops over both sites)
source("run_models/run_management_effect/run_management_effect_RF.R")   # + MLP, XGBoost, PI

# 3 — Figures            metrics_graphs.R first: it writes the CSVs the next two read
source("figures/metrics_graphs.R")
source("figures/metrics_for_individual_gaps_graphs.R")
source("figures/management_effect_graphs.R")
source("figures/VI_graphs.R")                          # independent — any order
source("figures/real_gaps_cleaveland_graphs.R")        # independent
source("figures/timeseries_nee_graphs.R")              # independent
source("figures/wind_rose_graphs.R")                   # independent
source("figures/montly_temperature_rainfall_table.R")  # independent — prints a table
```

---

## Folders

**`data_preparation/`** — quality control and predictor construction. `quality_control_nee/` screens each site-year's raw fluxes (range, instrument diagnostics, u\* threshold, stationarity flag, footprint), with the Kljun (2015) footprint model and boundary polygon in `grid_qc/`. Numbered scripts `01–09` then build the modelling dataset in order: `01` is the master runner that sources `02–07` for every site-year; `08` and `09` add the cross-validation gaps and the phytomass index. `ManagementEvents.R` defines the grazing/fertiliser tables and `generate_management_csvs.R` exports them (also called by `01`). `Regrowth_period.R` is an optional standalone annotation utility and is not part of the main pipeline.

**`models/`** — one cross-validation gap-filler per model. These are **not run directly**; a `run_models/` script loads the data and sources the matching `{MODEL}_CV.R`. Each trains once per artificial gap and predicts the withheld block. `miniRECgap_CV.R` produces Reco and GPP directly (NEE = Reco − GPP); `MDS_CV.R` fills NEE only and is excluded from L/VL gaps.

**`run_models/`** — entry points that load `{SITE}_cv.rds` and call a model. `run_ML_models/` holds the 12 machine-learning runs (site × model × managed/unmanaged); `run_benchmark_models/` holds MDS and miniRECgap; `run_management_effect/` runs the feature-addition study (BASE predictors + one management variable, looping over both sites).

**`figures/`** — every paper related figures and tables, built from the saved predictions and the prepared / QC data. `metrics_graphs.R` computes MAE/RMSE/R² (pooled and split by ≤30 d / >30 d since grazing) and writes `graphs/metrics_csv/`, which `metrics_for_individual_gaps_graphs.R` and `management_effect_graphs.R` consume; the remaining five scripts (`VI_graphs.R`, `real_gaps_cleaveland_graphs.R`, `timeseries_nee_graphs.R`, `wind_rose_graphs.R`, `montly_temperature_rainfall_table.R`) are independent and run in any order.

---

## Data availability

The eddy-covariance flux data, management records, and field boundary polygons (`nasco_site_list_wkt.csv`) are **not included** in this repository. They are available on reasonable request from **Ryan Burger** (ryan.burger@teagasc.ie).

---

## Data layout

Data files are **not included** (see [Data availability](#data-availability)). Place them under `data/`:

```text
data/                                          
├── raw_data/
│   ├── {SITE}_{2023,2024}_{eddypro,biomet,meta,fluxnet}.xlsx
│   └── JC1_2020_{eddypro,biomet,extra,AGC}.xlsx          
├── meteireann_data/
│   ├── met_hourly.csv                          hourly Temp, RH, rain
│   └── solar_hourly.xlsx                        hourly global radiation
├── canopy_data/
│   └── grass_drybiomass_JC_2023_2024.rds        date, height_{SITE}, drymass_{SITE}
├── management_data/                             written by generate_management_csvs.R
│   ├── grazing_events_{SITE}.csv
│   ├── fertiliser_events_{SITE}.csv
│   ├── last_events_before_year_{SITE}.csv
│   └── nitrogen_amounts_{SITE}.csv
├── data_qc/                                     Stage 0 output
│   └── results_qc_{SITE}_{YEAR}/
│       ├── merged_qc.{csv,rds}
│       ├── merged_qc_fullgrid.csv
│       └── ustar_seasonal_thresholds.rds
└── data_prepared/                               Stage 1 output
    ├── {SITE}_{YEAR}_predictors.rds
    ├── {SITE}.rds
    └── {SITE}_cv.rds                            ← model input
```

---

## Requirements

The scripts work with R version **R 4.6.1**.

```r
install.packages(c(
  "here", "tidyverse", "lubridate", "zoo", "hms", "readxl",
  "glue", "rlang", "scales",              # utilities (bundled with tidyverse)
  "patchwork", "ggh4x",                   # figure layout + faceting extensions
  "sf", "sp", "raster", "stars",          # spatial grid / footprint handling
  "ranger", "xgboost",                    # RF, XGBoost
  "REddyProc",                            # MDS benchmark + u* tools
  "reticulate", "tensorflow", "keras3"    # MLP (Python/TensorFlow backend)
))
```

The MLP scripts additionally need a working Python (v.3.10.12), TensorFlow (v.2.19.1) and Keras (v.3.10.0) environments behind `reticulate`.