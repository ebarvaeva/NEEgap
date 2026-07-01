# calc_footprint_FFP_mod.R — 2-D flux footprint model (Kljun et al. 2015)
#
# Defines calc_footprint_FFP_mod(): an R port of the Kljun FFP Python code. From the
# measurement height, roughness/wind, PBLH, Obukhov length, sigma_v, u* and wind
# direction it builds the crosswind-integrated footprint, spreads it laterally with a
# Gaussian, rotates it into the wind direction and shifts it onto the tower
# coordinates. Returns the density grid f_2d with its x_2d / y_2d coordinates.
#
# Input  : scalar micrometeorology for one half-hour (zm, z0 or umean, h, ol, sigmav, ustar, wind_dir)
# Output : list(x_ci_max, x_ci, f_ci, x_2d, y_2d, f_2d, flag_err). Sourced by the QC scripts; saves nothing.

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


  # --- Model parameters ---
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

  # --- Scaled parameters ---
  # Generate a sequence of scaled x* values (nx+2 values, then drop the first)
  xstar_ci_param <- seq(d, xstar_end, length.out = nx + 2)[-1]

  # Compute the scaled crosswind integrated footprint f*
  fstar_ci_param <- a * (xstar_ci_param - d)^b * exp(-c_param / (xstar_ci_param - d))

  # Scaled lateral dispersion parameter
  sigystar_param <- ac * sqrt(bc * xstar_ci_param^2 / (1 + cc * xstar_ci_param))

  # --- Determine the branch: use z0 if provided, otherwise umean ---
  if (!is.null(z0)) {
    # Compute psi_f based on the Obukhov length conditions
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

  # --- Compute the real-scale lateral dispersion (sigy) ---
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

  # --- Build the two-dimensional footprint ---
  # Assume uniform spacing in the x-direction; compute dx from x_ci
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

  # --- Apply rotation if a wind direction is specified ---
  if (!is.null(wind_dir)) {
    wind_rad <- wind_dir * pi / 180
    dist_mat <- sqrt(x_2d^2 + y_2d^2)
    angle_mat <- atan2(y_2d, x_2d)
    x_2d <- dist_mat * sin(wind_rad - angle_mat)
    y_2d <- dist_mat * cos(wind_rad - angle_mat)
  }

  # --- Shift the grid to the tower coordinates if provided ---
  if (!is.null(tower_x) && !is.null(tower_y)) {
    x_2d <- x_2d + tower_x
    y_2d <- y_2d + tower_y
  }

  # --- Optionally plot the footprint ---
  if (fig) {
    image(x_2d, y_2d, f_2d, main = "Flux Footprint", xlab = "x [m]", ylab = "y [m]", col = heat.colors(100))
    contour(x_2d, y_2d, f_2d, add = TRUE, col = "white")
  }

  # --- Return the results ---
  flag_err <- 0  # No error flag set

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
