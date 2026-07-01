# 08__artificial_gaps.R — construct artificial cross-validation gaps
#
# For each site, loads the prepared predictors, keeps rows with observed NEE,
# and builds four sets of contiguous leave-one-gap-out blocks (VL ~30 d,
# L ~14 d, M ~7 d, S ~1 d). For every gap label (e.g. M3) two columns are added:
# the logical flag (M3) and the masked training target (NEE_M3 = NEE_orig with
# gap rows set to NA). Gap structure matches the model CV scripts.
#
# Input  : data/data_prepared/{SITE}.rds
# Output : data/data_prepared/{SITE}_cv.rds

library(here)
library(dplyr)
library(tibble)
library(purrr)

SITES    <- c("JC1", "JC2")
IN_DIR   <- here::here("data", "data_prepared")
OUT_DIR  <- here::here("data", "data_prepared")

# target block size (in half-hours) and tolerance per gap-size class
gap_cfg <- list(
  VL = list(base = 30 * 48L, tol = 250L),
  L  = list(base = 14 * 48L, tol = 120L),
  M  = list(base =  7 * 48L, tol =  60L),
  S  = list(base =  1 * 48L, tol =  10L)
)

# split n rows into contiguous blocks ~base_size long (deviation <= tol)
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

# add flag columns (prefix1, prefix2, ...) and masked-NEE columns (NEE_prefix1, ...)
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

for (site in SITES) {
  
  # load prepared data (skip the site if missing)
  in_path <- file.path(IN_DIR, paste0(site, ".rds"))
  if (!file.exists(in_path)) {
    warning("Input not found, skipping: ", in_path); next
  }
  df_raw <- readRDS(in_path)
  
  # keep only rows with observed NEE
  df_cv <- df_raw %>% dplyr::filter(!is.na(NEE_orig))
  
  # build all four gap-size classes so their columns coexist in one data frame
  cv_VL <- create_cv_gaps(df_cv,        gap_cfg$VL$base, gap_cfg$VL$tol, "VL")
  cv_L  <- create_cv_gaps(cv_VL$data,   gap_cfg$L$base,  gap_cfg$L$tol,  "L")
  cv_M  <- create_cv_gaps(cv_L$data,    gap_cfg$M$base,  gap_cfg$M$tol,  "M")
  cv_S  <- create_cv_gaps(cv_M$data,    gap_cfg$S$base,  gap_cfg$S$tol,  "S")
  
  df_cv <- cv_S$data
  
  # save the CV data frame for this site
  out_path <- file.path(OUT_DIR, paste0(site, "_cv.rds"))
  saveRDS(df_cv, out_path)
}
