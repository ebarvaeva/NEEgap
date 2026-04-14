# =============================================================================
# SCRIPT 08 — Artificial Gap Construction
# =============================================================================
# PURPOSE
#   For each site (JC1, JC2), load the prepared predictor RDS, filter to rows
#   where NEE_orig is observed (non-NA), and construct four sets of contiguous
#   artificial gaps (VL / L / M / S) using leave-one-gap-out blocking.
#
#   The gap structure exactly mirrors the CV scheme used in RF_CV.R,
#   XGBoost_CV.R, and MLP_CV.R:
#
#     VL — Very Long : ~30 days  (30 × 48 half-hours), tolerance ±250 rows
#     L  — Long      : ~14 days  (14 × 48 half-hours), tolerance ±120 rows
#     M  — Medium    : ~ 7 days  ( 7 × 48 half-hours), tolerance ± 60 rows
#     S  — Short     : ~ 1 days  ( 1 × 48 half-hours), tolerance ± 10 rows
#
#   For each gap label (e.g. M3) two columns are added:
#     M3      — logical flag; TRUE where this gap falls
#     NEE_M3  — NEE_orig with M3 rows set to NA  (training target column)
#
# INPUTS
#   data/data_prepared/JC1.rds
#   data/data_prepared/JC2.rds
#
# OUTPUTS
#   data/data_prepared/JC1_cv.rds
#   data/data_prepared/JC2_cv.rds
#
#   Each output is the filtered input data frame augmented with:
#     VL1 … VLk   — gap flag columns (logical)
#     L1  … Lk
#     M1  … Mk
#     S1  … Sk
#     NEE_VL1 … NEE_VLk  — masked NEE columns (real, NA inside gap)
#     NEE_L1  … NEE_Lk
#     NEE_M1  … NEE_Mk
#     NEE_S1  … NEE_Sk
# =============================================================================

library(here)
library(dplyr)
library(tibble)
library(purrr)

# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------
SITES    <- c("JC1", "JC2")
IN_DIR   <- here::here("data", "data_prepared")
OUT_DIR  <- here::here("data", "data_prepared")

gap_cfg <- list(
  VL = list(base = 30 * 48L, tol = 250L),
  L  = list(base = 14 * 48L, tol = 120L),
  M  = list(base =  7 * 48L, tol =  60L),
  S  = list(base =  1 * 48L, tol =  10L)
)

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS  (identical to RF_CV.R)
# -----------------------------------------------------------------------------

# Divide n rows into contiguous blocks of approximately base_size rows.
# Remainder rows are distributed so that no block deviates from base_size
# by more than tol rows.
partition_timeseries_into_blocks <- function(n, base_size, tol) {
  stopifnot(n > 0, base_size >= 1, tol >= 0)
  k     <- max(1L, round(n / base_size))
  sizes <- rep(base_size, k)
  delta <- n - sum(sizes)
  if (delta != 0) {
    m    <- min(k, max(1L, ceiling(abs(delta) / max(1, tol))))
    step <- sign(delta) * floor(abs(delta) / m)
    rem  <- sign(delta) * (abs(delta) - abs(step) * m)
    sizes[(k - m + 1):k] <- sizes[(k - m + 1):k] + step
    if (rem != 0)
      sizes[(k - abs(rem) + 1):k] <- sizes[(k - abs(rem) + 1):k] + sign(rem)
  }
  for (i in k:2) {
    dev <- sizes[i] - base_size
    if (abs(dev) > tol) {
      ex           <- sign(dev) * (abs(dev) - tol)
      sizes[i]     <- base_size + sign(dev) * tol
      sizes[i - 1] <- sizes[i - 1] + ex
    }
  }
  sizes[sizes < 1] <- 1L
  sizes[k] <- sizes[k] + (n - sum(sizes))
  starts <- c(1L, head(cumsum(sizes) + 1L, -1L))
  ends   <- cumsum(sizes)
  stopifnot(ends[length(ends)] == n)
  list(sizes = sizes, indices = purrr::map2(starts, ends, seq.int))
}

# Build flag matrix (prefix1, prefix2, …) and masked-NEE matrix (NEE_prefix1, …).
# Returns a list with $data (augmented data frame) and $sizes (block sizes).
create_cv_gaps <- function(flux_data, base_size, tol, prefix) {
  n  <- nrow(flux_data)
  bl <- partition_timeseries_into_blocks(n,
                                         as.integer(round(base_size)),
                                         as.integer(round(tol)))
  fn <- paste0(prefix, seq_along(bl$indices))   # flag column names
  nn <- paste0("NEE_", fn)                       # masked-NEE column names
  
  fm <- matrix(FALSE,    nrow = n, ncol = length(fn), dimnames = list(NULL, fn))
  nm <- matrix(NA_real_, nrow = n, ncol = length(fn), dimnames = list(NULL, nn))
  
  for (i in seq_along(bl$indices)) {
    fm[bl$indices[[i]], i] <- TRUE
    nm[, i]                <- flux_data$NEE_orig
    nm[bl$indices[[i]], i] <- NA_real_
  }
  
  list(
    data  = dplyr::bind_cols(flux_data,
                             tibble::as_tibble(as.data.frame(fm)),
                             tibble::as_tibble(as.data.frame(nm))),
    sizes = bl$sizes
  )
}

# -----------------------------------------------------------------------------
# MAIN LOOP — process each site
# -----------------------------------------------------------------------------
for (site in SITES) {
  
  message("\n", strrep("=", 70))
  message("Processing site: ", site)
  message(strrep("=", 70))
  
  # --- 1. Load prepared data -------------------------------------------------
  in_path <- file.path(IN_DIR, paste0(site, ".rds"))
  if (!file.exists(in_path)) {
    warning("Input not found, skipping: ", in_path); next
  }
  df_raw <- readRDS(in_path)
  message("Loaded: ", in_path, "  (", nrow(df_raw), " rows)")
  
  # --- 2. Filter to rows with observed NEE -----------------------------------
  df_cv <- df_raw %>% dplyr::filter(!is.na(NEE_orig))
  message("Rows retained (NEE_orig not NA): ", nrow(df_cv),
          "  (dropped: ", nrow(df_raw) - nrow(df_cv), ")")
  
  # --- 3. Build all four gap-size categories ---------------------------------
  # Categories are stacked sequentially so that all gap columns coexist in
  # a single data frame, matching the structure expected by the model scripts.
  message("Constructing artificial gaps (VL → L → M → S) ...")
  
  cv_VL <- create_cv_gaps(df_cv,        gap_cfg$VL$base, gap_cfg$VL$tol, "VL")
  cv_L  <- create_cv_gaps(cv_VL$data,   gap_cfg$L$base,  gap_cfg$L$tol,  "L")
  cv_M  <- create_cv_gaps(cv_L$data,    gap_cfg$M$base,  gap_cfg$M$tol,  "M")
  cv_S  <- create_cv_gaps(cv_M$data,    gap_cfg$S$base,  gap_cfg$S$tol,  "S")
  
  df_cv <- cv_S$data
  
  # --- 4. Report gap counts --------------------------------------------------
  message("Gap counts:")
  message("  VL gaps : ", length(cv_VL$sizes),
          "  (sizes ", min(cv_VL$sizes), "–", max(cv_VL$sizes), " rows)")
  message("  L  gaps : ", length(cv_L$sizes),
          "  (sizes ", min(cv_L$sizes),  "–", max(cv_L$sizes),  " rows)")
  message("  M  gaps : ", length(cv_M$sizes),
          "  (sizes ", min(cv_M$sizes),  "–", max(cv_M$sizes),  " rows)")
  message("  S  gaps : ", length(cv_S$sizes),
          "  (sizes ", min(cv_S$sizes),  "–", max(cv_S$sizes),  " rows)")
  
  cat("VL gaps:", paste0("VL", seq_along(cv_VL$sizes), collapse = ", "), "\n")
  cat("L  gaps:", paste0("L",  seq_along(cv_L$sizes),  collapse = ", "), "\n")
  cat("M  gaps:", paste0("M",  seq_along(cv_M$sizes),  collapse = ", "), "\n")
  cat("S  gaps:", paste0("S",  seq_along(cv_S$sizes),  collapse = ", "), "\n")
  
  # --- 5. Report new columns added -------------------------------------------
  new_flag_cols <- c(
    paste0("VL", seq_along(cv_VL$sizes)),
    paste0("L",  seq_along(cv_L$sizes)),
    paste0("M",  seq_along(cv_M$sizes)),
    paste0("S",  seq_along(cv_S$sizes))
  )
  new_nee_cols <- paste0("NEE_", new_flag_cols)
  
  message("New columns added: ",
          length(new_flag_cols), " flag + ",
          length(new_nee_cols),  " masked-NEE = ",
          length(new_flag_cols) + length(new_nee_cols), " total")
  message("Output dimensions: ", nrow(df_cv), " rows × ", ncol(df_cv), " cols")
  
  # --- 6. Save ---------------------------------------------------------------
  out_path <- file.path(OUT_DIR, paste0(site, "_cv.rds"))
  saveRDS(df_cv, out_path)
  message("Saved: ", out_path)
}

message("\nScript 08 complete.")