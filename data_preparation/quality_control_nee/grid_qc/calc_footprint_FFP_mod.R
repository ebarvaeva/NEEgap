# =============================================================================
# calc_footprint_FFP_mod.R  -- Kljun et al. (2015) FFP Footprint Model (R port)
# =============================================================================
#
# PURPOSE
#   Computes a two-dimensional flux footprint probability density using the
#   Flux Footprint Prediction (FFP) parameterisation of Kljun et al. (2015).
#   Returns a 2D grid of footprint density values (f_2d) in real-world
#   coordinates centred at the tower and optionally rotated to the wind
#   direction.  The footprint grid is used by quality_control_nee.R to
#   compute the fraction of flux mass inside the site boundary polygon.
#
# MODEL OVERVIEW
#   The FFP model parameterises the crosswind-integrated footprint f*(x*) as
#   a scaled analytical function of the dimensionless distance x*:
#     f*(x*) = a * (x* - d)^b * exp(-c / (x* - d))
#   with empirical parameters a=1.4524, b=-1.9914, c=1.4622, d=0.1359.
#   The lateral crosswind distribution is Gaussian with dispersion sigma_y.
#   Physical (metre-scale) distances are recovered by scaling x* by zm and
#   a stability-dependent wind speed ratio.  Two branches are available:
#     z0 branch  -- uses roughness length and Obukhov stability correction
#     umean branch -- uses mean wind speed and u*/kappa
#
# COORDINATE SYSTEM
#   Footprint grid is generated in a tower-centred, flow-aligned frame
#   (x = downwind, y = crosswind).  If wind_dir is supplied, the grid is
#   rotated to align with geographic north so that points are in the same
#   projected CRS as the site boundary polygon (EPSG:2157 in this pipeline).
#   If tower_x and tower_y are given, the grid is also translated to absolute
#   map coordinates.
#
# INPUTS (all scalars for one half-hour)
#   zm       -- effective measurement height z - d  (m)
#   z0       -- roughness length (m); if NULL, umean must be provided
#   umean    -- mean wind speed (m/s); alternative to z0
#   h        -- planetary boundary layer height (m); from pblh_calc.R
#   ol       -- Obukhov stability length (m)
#   sigmav   -- std dev of lateral wind component (m/s)
#   ustar    -- friction velocity (m/s)
#   wind_dir -- wind direction (degrees from North); optional, used for rotation
#   nx       -- number of downwind grid points (default 200; increase for finer grids)
#   tower_x, tower_y -- projected tower coordinates (m) for absolute positioning
#
# RETURNS  (named list)
#   x_ci_max -- downwind distance of maximum crosswind-integrated footprint
#   x_ci, f_ci -- downwind distances and crosswind-integrated footprint values
#   x_2d, y_2d -- coordinate grids (m), shape (nx, n_y_total)
#   f_2d     -- 2D footprint probability density, same shape as x_2d/y_2d
#   flag_err -- 0 = no error
#
# REFERENCE
#   Kljun, N., Calanca, P., Rotach, M.W., Schmid, H.P. (2015). A simple
#   two-dimensional parameterisation for Flux Footprint Prediction (FFP).
#   Geoscientific Model Development, 8, 3695-3713.
#   https://doi.org/10.5194/gmd-8-3695-2015
#
# =============================================================================
calc_footprint_FFP_mod <- function(zm, z0 = NULL, umean = NULL, h, ol, sigmav, ustar,
                                   wind_dir = NULL,
                                   rs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
                                   rslayer = 0, nx = 200, crop = FALSE, fig = FALSE,
                                   tower_x = NULL, tower_y = NULL) {
  # Input validation
  if (is.null(zm) || is.null(h) || is.null(ol) || is.null(sigmav) || is.null(ustar) ||
      (is.null(z0) && is.null(umean))) {
    stop("FFP calculation: one or more required parameters are missing.")
  }


  # --- Model parameters (Kljun et al., 2015, Table 1) ----------------------
  # These are empirical constants fitted to LES simulations across a wide
  # range of stability conditions and surface roughnesses.
  a <- 1.4524
  b <- -1.9914
  c_param <- 1.4622  # named differently from R's reserved 'c'
  d <- 0.1359
  ac <- 2.17
  bc <- 1.66
  cc <- 20.0
  xstar_end <- 30
  oln <- 5000
  k <- 0.4

  # --- Dimensionless (scaled) footprint parameters --------------------------
  # x* is the dimensionless downwind distance; f*(x*) is the crosswind-
  # integrated footprint in scaled coordinates.
  # xstar_end = 30 defines the far-field cutoff (~99% of total footprint mass).
  xstar_ci_param <- seq(d, xstar_end, length.out = nx + 2)[-1]

  # Compute the scaled crosswind integrated footprint f*
  fstar_ci_param <- a * (xstar_ci_param - d)^b * exp(-c_param / (xstar_ci_param - d))

  # Scaled lateral dispersion parameter
  sigystar_param <- ac * sqrt(bc * xstar_ci_param^2 / (1 + cc * xstar_ci_param))

  # --- Convert from scaled to physical coordinates -------------------------
  # Two branches depending on whether roughness length (z0) or mean wind
  # speed (umean) is available.  The z0 branch uses the log-law stability
  # correction psi_f; the umean branch uses the measured wind profile directly.
  if (!is.null(z0)) {
    # Compute psi_f based on the Obukhov length conditions
    # Stability correction function psi_f (Paulson, 1970 / Brutsaert, 1982):
    #   Unstable/neutral (ol <= 0 or |ol| >= oln): Businger-Dyer correction
    #   Stable (ol > 0): linear approximation -5.3 * zm/ol
    if (ol <= 0 || ol >= oln) {
      xx <- (1 - 19.0 * zm / ol)^(0.25)
      psi_f <- log((1 + xx^2) / 2) + 2 * log((1 + xx) / 2) - 2 * atan(xx) + pi/2
    } else {
      psi_f <- -5.3 * zm / ol
    }
    x_ci <- xstar_ci_param * zm / (1 - (zm / h)) * (log(zm / z0) - psi_f)
    f_ci <- fstar_ci_param / zm * (1 - (zm / h)) / (log(zm / z0) - psi_f)
    xstarmax <- -c_param / b + d
    x_ci_max <- xstarmax * zm / (1 - (zm / h)) * (log(zm / z0) - psi_f)
  } else {
    x_ci <- xstar_ci_param * zm / (1 - (zm / h)) * (umean / ustar * k)
    f_ci <- fstar_ci_param / zm * (1 - (zm / h)) / (umean / ustar * k)
    xstarmax <- -c_param / b + d
    x_ci_max <- xstarmax * zm / (1 - (zm / h)) * (umean / ustar * k)
  }

  # --- Real-scale lateral dispersion sigma_y --------------------------------
  # scale_const adjusts the lateral spread based on stability; it accounts for
  # the increase in lateral diffusion under unstable conditions relative to
  # neutral.  Values are capped at 1.0 to prevent over-correction.
  if (abs(ol) > oln) {
    ol <- -1e6
  }
  if (ol <= 0) {
    scale_const <- 1e-5 * (1 / abs(zm / ol)) + 0.80
  } else {
    scale_const <- 1e-5 * (1 / abs(zm / ol)) + 0.55
  }
  if (scale_const > 1) {
    scale_const <- 1.0
  }
  sigy <- sigystar_param / scale_const * zm * sigmav / ustar
  sigy[sigy < 0] <- NA

  # --- Construct the 2D footprint grid -------------------------------------
  # For each downwind position x_ci[i], the crosswind profile is a Gaussian
  # with standard deviation sigy[i].  f_pos holds the positive-y half;
  # f_neg is its mirror image.  They are combined into f_2d.
  if (length(x_ci) >= 2) {
    dx <- x_ci[2] - x_ci[1]
  } else {
    dx <- 1
  }
  # Create y positions for the positive half (from 0 upward)
  y_pos <- seq(0, (length(x_ci) / 2) * dx * 1.5, by = dx)

  n_rows <- length(f_ci)
  n_cols <- length(y_pos)
  f_pos <- matrix(NA, nrow = n_rows, ncol = n_cols)

  for (i in 1:n_rows) {
    # For each row, compute the crosswind profile assuming a Gaussian shape
    f_pos[i, ] <- f_ci[i] * (1 / (sqrt(2 * pi) * sigy[i])) *
      exp(- (y_pos^2) / (2 * sigy[i]^2))
  }

  # Mirror to get negative y-values
  y_neg <- -rev(y_pos)
  # To avoid duplication at zero, remove the duplicate (i.e. remove the first element)
  y_neg <- y_neg[-1]

  y_combined <- c(y_neg, y_pos)

  # Mirror f_pos horizontally to obtain the negative side (drop the duplicate column)
  f_neg <- f_pos[, ncol(f_pos):1]
  f_neg <- f_neg[, -ncol(f_neg), drop = FALSE]

  # Concatenate to form the full 2D footprint grid
  f_2d <- cbind(f_neg, f_pos)

  # Create corresponding x and y grids.
  n_total <- length(y_combined)
  # x_2d: replicate the x_ci vector across columns
  x_2d <- matrix(rep(x_ci, each = n_total), nrow = n_rows, ncol = n_total, byrow = TRUE)
  # y_2d: replicate the combined y vector across rows
  y_2d <- matrix(rep(y_combined, times = n_rows), nrow = n_rows, ncol = n_total, byrow = TRUE)

  # --- Rotate from flow-aligned to geographic coordinates ------------------
  # wind_dir (degrees from North, clockwise) specifies the direction FROM which
  # the wind is blowing.  Rotating by wind_dir aligns the footprint with the
  # geographic coordinate system so that x_2d/y_2d are in the same CRS as
  # the site boundary polygon (EPSG:2157).
  if (!is.null(wind_dir)) {
    wind_rad <- wind_dir * pi / 180
    dist_mat <- sqrt(x_2d^2 + y_2d^2)
    angle_mat <- atan2(y_2d, x_2d)
    x_2d <- dist_mat * sin(wind_rad - angle_mat)
    y_2d <- dist_mat * cos(wind_rad - angle_mat)
  }

  # --- Translate to absolute map coordinates (tower-relative -> world) -----
  if (!is.null(tower_x) && !is.null(tower_y)) {
    x_2d <- x_2d + tower_x
    y_2d <- y_2d + tower_y
  }

  # --- Optionally plot the footprint ---
  if (fig) {
    image(x_2d, y_2d, f_2d, main = "Flux Footprint", xlab = "x [m]", ylab = "y [m]", col = heat.colors(100))
    contour(x_2d, y_2d, f_2d, add = TRUE, col = "white")
  }

  # --- Return results -------------------------------------------------------
  flag_err <- 0  # 0 = successful computation; non-zero values reserved for future errors

  return(list(
    x_ci_max = x_ci_max,
    x_ci = x_ci,
    f_ci = f_ci,
    x_2d = x_2d,
    y_2d = y_2d,
    f_2d = f_2d,
    flag_err = flag_err
    # Note: For brevity, contour level (rs) outputs are not computed.
  ))
}
