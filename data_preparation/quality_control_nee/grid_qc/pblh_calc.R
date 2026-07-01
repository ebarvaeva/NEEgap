# pblh_calc.R — planetary boundary layer height per half-hour (Kljun 2015)
#
# Defines boundary_layer_height(): from Obukhov length, friction velocity, sigma_v,
# latitude, (z-d)/L and air temperature it returns the PBLH used by the 2-D footprint
# model. Stable half-hours ((z-d)/L >= 0) use the diagnostic equilibrium height;
# unstable half-hours integrate the growth-rate equation forward from the last stable
# value. Any half-hour with a missing input returns NA.
#
# Input  : numeric vectors, one value per half-hour (Ls, ustars, t_covs, LAT, zLs, air_ts)
# Output : numeric vector of PBLH in metres. Sourced by the QC scripts; saves nothing.

boundary_layer_height <- function(Ls, ustars, t_covs, LAT, zLs, air_ts) {
  # Constants
  omg <- 7.2921159e-5
  # Calculate the Coriolis parameter (f) using latitude in degrees.
  f <- 2.0 * omg * sin(LAT * pi/180)
  neutral_limit <- 0
  A <- 0.2
  B <- 2.5
  C <- 8
  gamma <- 0.01
  g <- 9.81
  
  # Initialize the output vector.
  hs <- rep(NA, length(Ls))
  
  # Calculate the stable boundary layer height component
  h_stab <- (Ls / 3.8) * (-1.0 + sqrt(1.0 + 2.28 * (ustars / (f * Ls))))
  
  first_valid_idx <- which(!is.na(h_stab))[1]
  h_cur <- NA
  
  for (i in seq_along(Ls)) {
    L_val <- Ls[i]
    zL_val <- zLs[i]
    ustar_val <- ustars[i]
    t_cov_val <- t_covs[i]
    air_t_val <- air_ts[i]
    
    # If any of the required inputs are missing, assign NA and continue.
    if (any(is.na(c(L_val, zL_val, ustar_val, t_cov_val, air_t_val)))) {
      hs[i] <- NA
      next
    }
    
    if (zL_val >= neutral_limit) {
      # Stable conditions: use the stable formulation.
      h_calc <- (L_val / 3.8) * (-1.0 + sqrt(1.0 + 2.28 * (ustar_val / (f * L_val))))
      hs[i] <- h_calc
      h_cur <- h_calc
      next
    } else {
      # Unstable conditions: adjust the current boundary layer height.
      if (is.na(h_cur)) {
        if (!is.na(first_valid_idx)) {
          h_cur <- h_stab[first_valid_idx]
        }
        if (is.na(h_cur)) {
          h_cur <- 100.0
        }
      }
      D <- t_cov_val / gamma
      denom1 <- ((1 + 2 * A) * h_cur) - (2 * B * 0.4 * L_val)
      if (denom1 <= 0) {
        hs[i] <- NA
        next
      }
      E <- (h_cur^2) / denom1
      denom2 <- ((1 + A) * h_cur) - (B * 0.4 * L_val)
      if (denom2 <= 0) {
        hs[i] <- NA
        next
      }
      F <- (C * (ustar_val^2) * air_t_val) / (gamma * g * denom2)
      dh_dt <- (D / (E + F)) * 1800.0
      new_h_cur <- h_cur + dh_dt
      if (new_h_cur < 0) {
        new_h_cur <- 0
      }
      h_cur <- new_h_cur
      hs[i] <- h_cur
    }
  }
  return(hs)
}
