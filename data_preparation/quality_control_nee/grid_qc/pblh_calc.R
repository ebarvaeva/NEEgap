# =============================================================================
# pblh_calc.R  — Planetary Boundary Layer Height (PBLH)
# =============================================================================
#
# PURPOSE
#   Computes the planetary boundary layer height (PBLH, m) for each half-hour
#   using the hybrid stable/unstable formulation described in Kljun et al.
#   (2015) and its Python reference implementation.  PBLH is required as an
#   upper boundary condition by the FFP flux footprint model.
#
# ALGORITHM
#   Two stability regimes are handled:
#
#   STABLE (zL >= 0, i.e. (z-d)/L >= 0):
#     A diagnostic formula based on the Obukhov length and friction velocity:
#       h = (L / 3.8) * (-1 + sqrt(1 + 2.28 * u* / (f * L)))
#     where f = 2 * omega * sin(lat) is the Coriolis parameter.
#     This is a steady-state solution; each half-hour is independent.
#
#   UNSTABLE (zL < 0):
#     A prognostic (time-stepping) formulation.  The boundary layer grows
#     by dh/dt proportional to the convective velocity scale (sigma_v) and
#     the buoyancy flux, resisted by the stable temperature gradient above:
#       dh/dt = D / (E + F) * 1800
#     where D = sigma_v / gamma, E and F depend on h_cur, L, and u*.
#     h_cur carries state forward from the previous half-hour; the first
#     unstable half-hour initialises from the most recent stable PBLH.
#     Because of this sequential dependence, the function must be called
#     in row order (which quality_control_nee.R ensures via a for-loop).
#
# INPUTS
#   Ls     -- Obukhov length (m), scalar per call
#   ustars -- friction velocity u* (m/s)
#   t_covs -- standard deviation of lateral wind (sigma_v, m/s)
#   LAT    -- site latitude (degrees); 52.0 is a reasonable default for Ireland
#   zLs    -- stability parameter (z-d)/L; determines stable vs unstable branch
#   air_ts -- air temperature (K); used in the buoyancy denominator F
#
# CONSTANTS
#   omega = 7.2921e-5 rad/s  (Earth rotation rate)
#   A=0.2, B=2.5, C=8        (Kljun et al. empirical coefficients)
#   gamma = 0.01 K/m         (free-atmosphere lapse rate)
#   g = 9.81 m/s2
#   kappa = 0.4              (von Karman constant, used implicitly via B*kappa)
#
# RETURNS
#   Numeric scalar (or vector if called with vectors): PBLH in metres.
#   NA is returned for any row where an input is missing or a denominator
#   is non-positive (numerical guard to prevent Inf/NaN propagation).
#
# REFERENCE
#   Kljun, N., Calanca, P., Rotach, M.W., Schmid, H.P. (2015). A simple two-
#   dimensional parameterisation for Flux Footprint Prediction (FFP).
#   Geoscientific Model Development, 8, 3695-3713.
#
# =============================================================================
boundary_layer_height <- function(Ls, ustars, t_covs, LAT, zLs, air_ts) {
  # Constants (Kljun et al., 2015)
  omg <- 7.2921159e-5         # Earth's angular velocity (rad/s)
  f <- 2.0 * omg * sin(LAT * pi/180)  # Coriolis parameter (rad/s)
  neutral_limit <- 0          # zL threshold; >= 0 -> stable branch
  A <- 0.2                    # entrainment coefficient
  B <- 2.5                    # empirical coefficient in stable-layer resistance
  C <- 8                      # buoyancy flux coefficient
  gamma <- 0.01               # free-atmosphere temperature lapse rate (K/m)
  g <- 9.81                   # gravitational acceleration (m/s2)
  
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
      # STABLE: diagnostic formula; each row is independent.
      h_calc <- (L_val / 3.8) * (-1.0 + sqrt(1.0 + 2.28 * (ustar_val / (f * L_val))))
      hs[i] <- h_calc
      h_cur <- h_calc
      next
    } else {
      # UNSTABLE: prognostic; h_cur carries forward from the previous half-hour.
      if (is.na(h_cur)) {
        if (!is.na(first_valid_idx)) {
          h_cur <- h_stab[first_valid_idx]
        }
        if (is.na(h_cur)) {
          h_cur <- 100.0
        }
      }
      D      <- t_cov_val / gamma      # convective velocity scale term
      denom1 <- ((1 + 2 * A) * h_cur) - (2 * B * 0.4 * L_val)  # E denominator
      if (denom1 <= 0) {
        hs[i] <- NA
        next
      }
      E      <- (h_cur^2) / denom1     # shape factor
      denom2 <- ((1 + A) * h_cur) - (B * 0.4 * L_val)  # F denominator
      if (denom2 <= 0) {
        hs[i] <- NA
        next
      }
      F     <- (C * (ustar_val^2) * air_t_val) / (gamma * g * denom2)  # buoyancy term
      dh_dt <- (D / (E + F)) * 1800.0  # growth rate * 1800 s per half-hour
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
