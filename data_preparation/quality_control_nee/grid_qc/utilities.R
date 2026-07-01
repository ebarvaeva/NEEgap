# utilities.R — spatial-grid and footprint helper functions
#
# A library of helpers for building and post-processing 2-D flux footprints:
# atmospheric-stability class, grid construction / resampling, matrix<->long
# reshaping, cumulative-weight source-area masking, and the site helpers that read a
# boundary polygon and compute the fraction of footprint mass falling inside it.
# Defines functions only and saves nothing; sourced wherever these helpers are used.
#
# Input  : none (function definitions)
# Output : none

atmospheric_stability <- function(mo_length, z, zd, zo, neutral_thr = 0.04) {
  zm <- z - zd
  if (!missing(zo)) {
    zu <- zm * (log(zm / zo) - 1 + (zo / zm))
    zl <- zu / mo_length
  } else {
    zl <- zm / mo_length
  }
  dplyr::case_when(
    abs(zl) < neutral_thr ~ "neutral",
    zl < 0 ~ "unstable",
    zl > 0 ~ "stable"
  )
}

#' Coerce a gridded object to a matrix
#'
#' @param x A gridded object. Currently methods exist for classes data.frame,
#'   Raster, and stars.
#' @param ... Additional arguments to pass to coercion function.
#'
#' @return A matrix.
#'
#' @importFrom methods "as"
#'
#' @export
as_matrix <- function(x, ...) {

  if (inherits(x, "data.frame")) {
    return(as.matrix(x, ...))
  }

  if (inherits(x, "stars")) {
    x <- as(x, "Raster")
  }

  raster::as.matrix(x, ...)
}

# Performs vector operations on matrices while retaining matrix dimensions
with_matrix <- function(.x, .f, ..., .dims = dim(.x)) {

  mat_fun <- rlang::as_function(.f)

  vec <- mat_fun(as.vector(.x), ...)

  matrix(vec, nrow = .dims[1], ncol = .dims[2])
}

with_matrix2 <- function(.x, .y, .f, ..., .dims = dim(.x)) {

  # if (!isTRUE(all.equal(.dims, dim(.y)))) {
  #   stop("Dimensions of .x and .y must be equal.", call. = FALSE)
  # }

  mat_fun <- rlang::as_function(.f)

  vec <- mat_fun(as.vector(.x), as.vector(.y), ...)

  matrix(vec, nrow = .dims[1], ncol = .dims[2])
}

#' Construct matrix template for calculating footprint
#'
#' @param extent A numeric vector of length four
#' @param fetch An integer, maximum length in each direction to calculate
#'   footprint probabilities. Standard is 100 * (z - zd).
#' @param res An integer, width of grid cells to be used in calculating
#'   footprint probabilities. Defaults to one (i.e. a 1-m resolution).
#'
#' @export
grid_init <- function(extent, fetch, res = 1) {

  dx <- res
  dy <- res

  if (missing(extent) & missing(fetch)) {
    stop("Must provide either extent limits or fetch.", call. = FALSE)
  } else if (!missing(extent)) {
    nx <- (extent[2] - extent[1]) / dx
    ny <- (extent[4] - extent[3]) / dy
  } else if (!missing(fetch)) {
    fetch <- floor(fetch)
    nx <- fetch * 2
    ny <- fetch * 2
    extent <- c(-nx * dx / 2, nx * dx / 2, -ny * dy / 2, ny * dy / 2)
  }

  # Put extent into more convenient vars
  xmin <- extent[1]
  xmax <- extent[2]
  ymin <- extent[3]
  ymax <- extent[4]

  # Define physical domain in cartesian coordinates
  x <- seq(xmin, xmax, length.out = nx + 1)
  y <- -seq(ymin, ymax, length.out = ny + 1)

  x_2d <- matrix(rep(x, each = ny + 1), nrow = ny + 1)
  y_2d <- matrix(rep(y, nx + 1), nrow = ny + 1)

  out <- list(x = x_2d, y = y_2d)
  attributes(out) <- list("names" = names(out), "fetch" = fetch, "res" = res)
  out
}

#' Resample a raster to a standard matrix grid
#'
#' @param x A Raster object to be resampled.
#' @param grid A list of matrices as returned by [grid_init()].
#' @param coords A numeric vector of length two indicating the spatial
#'   coordinates of the "center" grid cell.
#' @param method Name of the resampling method passed to raster::resample. Can
#'   be "ngb" (the default) or "bilinear".
#'
#' @return A matrix that matches spatial dimensions of grid.
#'
#' @importFrom sf "st_bbox"
#'
#' @export
snap_to_grid <- function(x, grid, coords, method = "ngb") {

  temp <- grid_template_rst(grid, coords)

  # Return the original object if extent is the same as grid
  # - tolerance is ~10 cm
  equal_bbox <- all.equal(
    as.numeric(sf::st_bbox(temp)), as.numeric(sf::st_bbox(x))
  )
  if (isTRUE(equal_bbox)) {
    return(as_matrix(x))
  }

  # Assign CRS from original object
  raster::crs(temp) <- raster::crs(x)

  # Resample according to grid parameters
  resampled <- raster::resample(x, temp, method = method)

  out <- as_matrix(resampled)

  # Mask AOI if exists in grid list (maybe not, easier if there are no NAs)
  if ("aoi" %in% names(grid)) {
    aoi_mask <- grid$aoi
    aoi_mask[which(aoi_mask == 0)] <- NA
    #out <- out * aoi_mask
  }

  # Return matrix object fit to grid
  out
}

aoi_to_grid <- function(aoi, coords, delta = 1, center = FALSE) {

  # Ensure that grid template will cover entire AOI
  # - square buffer of [length] around tower should contain AOI regardless of
  #   exact tower location
  bbox <- sf::st_bbox(aoi)
  length <- ceiling(max(bbox$xmax - bbox$xmin, bbox$ymax - bbox$ymin))

  # Construct oversized grid template
  template <- buffer_grid(coords, length, delta = delta, center = center)

  # Crop to actual AOI extent
  grid_crop <- template %>%
    sf::st_set_crs(sf::st_crs(aoi)) %>%
    sf::st_crop(bbox, as_points = FALSE) %>%
    tibble::as_tibble()

  # Set values as x coordinates (x grid)
  grid_x <- grid_crop %>%
    dplyr::mutate(values = x - coords[1]) %>%
    stars::st_as_stars() %>%
    as_matrix()

  # Set values as y coordinates (y grid)
  grid_y <- grid_crop %>%
    dplyr::mutate(values = y - coords[2]) %>%
    stars::st_as_stars() %>%
    as_matrix()

  # Return list of x and y grids
  list(x = grid_x, y = grid_y)
}

aoi_to_mask <- function(aoi, coords, delta = 1, center = FALSE) {

  # Ensure that grid template will cover entire AOI
  # - square buffer of [length] around tower should contain AOI regardless of
  #   exact tower location
  bbox <- sf::st_bbox(aoi)
  length <- ceiling(max(bbox$xmax - bbox$xmin, bbox$ymax - bbox$ymin))

  # Construct oversized grid template
  template <- buffer_grid(coords, length, delta = delta, center = center)

  # Crop to actual AOI extent
  grid_crop <- template %>%
    sf::st_set_crs(sf::st_crs(aoi)) %>%
    sf::st_crop(bbox, as_points = FALSE) %>%
    sf::st_crop(aoi)

  # Convert to matrix format
  grid_mat <- as_matrix(grid_crop)

  outside <- which(is.na(grid_mat))
  grid_mat[outside] <- 0

  grid_mat
}

cut_grid <- function(x, width, keep_partial = TRUE, coords = NULL) {

  # For creating labeled grid quadrants

  grid <- x
  if (inherits(x, "sf")) {
    if (is.null(coords)) {
      stop("'coords' must be provided if x is a spatial object.")
    }
    grid <- aoi_to_grid(x, coords)
  }

  #dim <- dim(grid[[1]])
  center <- width / 2

  # if (!is.null(center)) {
  #   length <- max(dim[1], dim[2]) %>% ceiling()
  #   template <- buffer_grid(center, length = length)
  # }

  parse_cut_label <- function(x, fun = mean) {

    num <- x %>%
      as.character() %>%
      stringr::str_remove_all("\\[|\\]|\\(|\\)") %>%
      stringr::str_split(",") %>%
      purrr::map(as.numeric)

    num %>%
      purrr::map(rlang::as_function(fun)) %>%
      purrr::simplify()
  }

  row_bins <- grid %>%
    purrr::pluck(2) %>%
    magrittr::extract(, 1) %>%
    ggplot2::cut_width(width = width, center = center) %>%
    parse_cut_label() %>%
    dplyr::if_else(. < 0, abs(. - 1), .) %>%
    as.factor() %>%
    as.integer()

  col_bins <- grid %>%
    purrr::pluck(1) %>%
    magrittr::extract(1, ) %>%
    ggplot2::cut_width(width = width, center = center) %>%
    parse_cut_label() %>%
    dplyr::if_else(. < 0, abs(. - 1), .) %>%
    as.factor() %>%
    as.integer()

  # Determine bin order based on distance from center
  bin_order <- expand.grid(row_bins, col_bins) %>%
    dplyr::distinct(Var1, Var2) %>%
    dplyr::mutate(
      code = paste0(
        stringr::str_pad(Var1, 2, pad = "0"),
        stringr::str_pad(Var2, 2, pad = "0")
      ),
      dist = sqrt(Var1^2 + Var2^2)
    ) %>%
    dplyr::arrange(dist) %>%
    dplyr::pull(code)

  #row_bins <- ggplot2::cut_width(grid[[2]][, 1], width = width, center = center)
  #col_bins <- ggplot2::cut_width(grid[[1]][1, ], width = width, center = center)

  # row_bins <- ggplot2::cut_number(seq(1, dim[1]), n, labels = FALSE)
  # col_bins <- ggplot2::cut_number(seq(1, dim[2]), n, labels = FALSE)
  grid_bins <- outer(
    stringr::str_pad(row_bins, 2, pad = "0"),
    stringr::str_pad(col_bins, 2, pad = "0"),
    FUN = paste0
  )
  grid_lab <-  with_matrix(
    grid_bins, ~ as.integer(forcats::fct_relevel(factor(.x), bin_order))
  )

  if (!keep_partial) {
    bin_cells <- table(grid_lab)
    partial <- which(bin_cells < width^2)
    grid_lab[grid_lab %in% partial] <- NA
  }

  grid_lab
}

get_trim_extent <- function(x, values = 0) {

  # Only support for one trim value at the moment
  keep_y <- which(rowSums(x) != values)
  keep_x <- which(colSums(x) != values)

  # ymin, ymax, xmin, xmax
  c(min(keep_y), max(keep_y), min(keep_x), max(keep_x))
}

trim_matrix <- function(x, extent = get_trim_extent(x)) {

  x[extent[1]:extent[2], extent[3]:extent[4]]
}

#' Pivot matrix data from wide to long
#'
#' @param x A matrix
#'
#' @importFrom tidyr "pivot_longer" tibble "as_tibble"
#'
#' @export
pivot_matrix <- function(x) {

  if (storage.mode(x) == "character") {
    x <- with_matrix(x, ~ unclass(as.factor(.x)))
  }

  dims <- dim(x)
  if (is.null(dimnames(x))) dimnames(x) <- list(rev(1:dims[1]), 1:dims[2])

  x %>%
    tibble::as_tibble(rownames = "y") %>%
    tidyr::pivot_longer(-y, names_to = "x", values_to = "z") %>%
    dplyr::mutate(dplyr::across(is.character, as.numeric))
}

#' Write a matrix to a delimited file
#'
#' @param x A matrix object.
#' @param path A path to write to.
#' @param trunc Number of decimal places used to integerize data contents. If
#'   set to 0 (the default), no changes are made. A common value is 9.
#' @param compress Use gzip compression? Must be either TRUE or FALSE. If TRUE
#'   (recommended), a ".gz" extension will be added to the end of file.
#'
#' @importFrom readr "write_delim"
#' @importFrom tibble "as_tibble"
#'
#' @export
write_matrix <- function(x, path, trunc = 0, compress = TRUE) {

  # Convert data type to integer for more efficient storage
  # - this means weights below (1 / 'trunc') are effectively zero
  if (trunc > 0) {
    x <- trunc(x * 10^trunc)
    storage.mode(x) <- "integer"
  }

  # Set path
  file_path <- paste0(path, ".txt")
  if (compress) file_path <- paste0(file_path, ".gz")

  # Convert to tbl for writing
  tbl <- tibble::as_tibble(x, .name_repair = "minimal")

  # Write to file
  readr::write_delim(tbl, file_path, col_names = FALSE)
}

#' Read whitespace-separated columns into a matrix
#'
#' @param file A path to a file.
#' @param trunc Number of decimal places used to integerize data contents. If
#'   contents were not truncated, set to 0 (the default). A common value is 9.
#'
#' @importFrom readr "read_table2" "cols" "col_integer" "col_double"
#'
#' @export
read_matrix <- function(file, trunc = 0) {

  if (trunc > 0) {
    data_type <- readr::col_integer()
  } else {
    data_type <- readr::col_double()
  }

  data <- readr::read_table2(
    file, col_names = FALSE, col_types = readr::cols(.default = data_type),
    progress = FALSE
  )
  # data <- vroom::vroom(
  #   file, delim = " ", col_names = FALSE,
  #   col_types = readr::cols(.default = data_type), progress = FALSE
  # )

  data_mat <- as_matrix(data)
  dimnames(data_mat) <- NULL

  # Scale based on how data was stored
  data_mat <- data_mat / 10^trunc

  data_mat
}

#' Read grid coordinates from two .txt files
#'
#' probably get rid of this - pretty unnecessary
#'
#' @param path A character vector containing a single path
#' @param names A character vector of length two, the file names of grid
#'   elements (omitting the extension)
#'
#' @importFrom magrittr "%>%"
#'
#' @export
read_grid <- function(path, names = c("x", "y")) {

  file_names <- paste0(names, ".txt")
  grid_files <- file.path(path, file_names)

  grid_files %>%
    purrr::map(read_matrix, trunc = 0) %>%
    rlang::set_names(names)
}

pluck_cell <- function(x, grid, coords) {

  cell_loc <- which(grid$x == coords[1] & grid$y == coords[2])
  x[cell_loc]
}

#' Rasterize a footprint matrix
#'
#' @param x A footprint object
#' @param grid A list of length two containing matrices of equal dimensions,
#'   indicating x and y coordinates. Template returned by \link{grid_init}.
#' @param coords A numeric vector of length two
#' @param crs A character or object of class CRS.
#'
#' @export
fp_rasterize <- function(x, grid, coords, crs) {

  # Add geographical coordinates to grid
  grid$y <- grid$y + coords[1] # lat
  grid$x <- grid$x + coords[2] # lon

  # Convert list of matrices to XYZ data frame
  fp_grid <- as.matrix(data.frame(
    # Necessary order: lat (y), lon (x), values
    y = as.vector(grid$y),
    x = as.vector(grid$x),
    z = as.vector(x)
  ))

  # Rasterize
  raster::rasterFromXYZ(fp_grid, crs = crs)

}

rotate_grid <- function(grid, dir) {
  fpd <- sqrt(grid$y^2 + grid$x^2)
  fpa <- atan2(grid$x, grid$y) * 180 / pi - dir

  out$x <- cos(fpa * pi/180) * fpd
  out$y <- -1 * sin(fpa * pi/180) * fpd
  out
}

accumulate_weights <- function(x, max = 1, zero_as_max = FALSE) {

  source <- x
  if (!zero_as_max) source <- with_matrix(x, ~ dplyr::na_if(.x, 0))

  # Sort matrix values from high to low
  sorted <- with_matrix(
    source, ~ vctrs::vec_sort(.x, direction = "desc", na_value = "smallest")
  )
  sort_order <- vctrs::vec_order(
    as.vector(source), direction = "desc", na_value = "smallest"
  )

  # Calculate cumulative sums of footprint weights
  summed <- source
  summed[sort_order] <- with_matrix(sorted, ~ cumsum(.x))

  accum_source <- with_matrix(summed, ~ tidyr::replace_na(.x, 0))

  accum_source[accum_source > max] <- 0

  accum_source
}

#' Mask footprint to its analytical source area
#'
#' As the analytical footprint approaches 100%, the coverage increases rapidly
#' over area that contributes little to the total measured flux (Matthes et al.,
#' 2014). Thus, it may be desirable to consider a percantage of the calculated
#' footprint as the analytical footprint. This is acheived by masking the
#' footprint to a given maximum cumulative weight.
#'
#' @param x The footprint as a numeric matrix or RasterLayer.
#' @param p A numeric value between 0 and 1, recommended by Kljun et al. (2015)
#'   to be less than 0.9. Defaults to 0.85.
#' @param mask_value A value to assign cells beyond the maximum percent weight.
#'   Defaults to NA, but 0 may be useful in some situations.
#'
#' @return An object of the same class as x, either a numeric matrix or
#'   RasterLayer.
#'
#' @export
mask_source_area <- function(x, p = 0.85, mask_value = NA) {

  # Assume that p is given as a percentage if >1
  if (p > 1) p <- p / 100

  # Sort matrix values from high to low
  x_vec <- as.vector(x)
  sorted <- vctrs::vec_sort(x_vec, direction = "desc", na_value = "smallest")

  # Calculate cumulative sums of footprint weights
  summed <- cumsum(sorted)

  # Get the minimum value of cumulative sum less than p
  source_area <- sorted[summed < p]

  # If cumulative footprint is never less than p, source is too close to tower
  if (length(source_area) == 0) {
    cutoff <- max(x)
  } else {
    cutoff <- min(source_area)
  }

  out <- x

  # Set all other cells to fill value
  out[out < cutoff] <- mask_value

  out
}

#' Convert footprint to a polygon representing a given source area
#'
#' @param x The footprint as a numeric matrix or RasterLayer.
#' @param p A numeric value between 0 and 1, recommended by Kljun et al. (2015)
#'   to be less than 0.9. Defaults to 0.85.
#'
#' @export
source_area_polygon <- function(x, p = 0.85) {

  masked <- mask_source_area(x, p = p, mask_value = NA)

  class <- class(x)

  if (class == "RasterLayer") {
    # Method for RasterLayer class footprint
    # Set all raster values to 0 - they don't matter anymore
    homog <- raster::calc(x, function(x) x * 0)
    # Return polygon of the total area
    raster::rasterToPolygons(homog, dissolve = TRUE)
  } else if (class == "SpatialLinesDataFrame") {
    # Method for SpatialLinesDataFrame class footprint
    # Tidy ftp for creating polygon object
    fort <- broom::tidy(x)
    # Create spatial polygons object
    poly <- sp::Polygon(fort[, 1:2])
    polys <- sp::Polygons(list(sp::Polygon(poly)), "poly")
    sp::SpatialPolygons(
      list(polys), proj4string = sp::CRS("+proj=utm +zone=18 +ellps=GRS80")
    )
  }
}

# ----------------------------------------------------------------------------------------------------------
# Function to extract footprint boundary and tower coordinates for a given site
# ----------------------------------------------------------------------------------------------------------
# Input:
#   - site: Site code ("BCN" or "CD")
#   - polygons_df: A dataframe containing site metadata, including boundary polygons and tower locations
#
# Output:
#   - A list containing:
#       - `site_name`: Short code of the site
#       - `full_name`: Full descriptive name of the site
#       - `tower_x`, `tower_y`: Tower coordinates (EPSG:2157)
#       - `boundary_poly`: Boundary polygon as an `sf` object
#       - `raw_polygon`: The original filtered row from `polygons_df`
#
# Notes:
#   - Site codes are mapped to descriptive names using an internal dictionary.
#   - Function validates site code and ensures polygon data is available.
#   - Coordinates are parsed from `tower_loc_2157` column (POINT geometry).

process_footprint_data <- function(site, polygons_df) {
  site_mapping <- list(
    BCN = "Ballycanew (ACP)",
    CD = "Castledockrell Arable (ACP)",
    JC1 = "Johnstown Castle - Teagasc (1)",
    JC2 = "Johnstown Castle - Teagasc (2)"
  )
  
  if (!site %in% names(site_mapping)) {
    stop("Invalid site code. Must be 'BCN' or 'CD'")
  }
  
  site_poly <- polygons_df |> 
    filter(`Site Name` == site_mapping[[site]])
  
  if (nrow(site_poly) == 0) {
    stop("Polygon data not found for site: ", site)
  }
  
  # Extract coordinates
  coords <- as.numeric(strsplit(
    gsub("[POINT()]", "", site_poly$tower_loc_2157), " ")[[1]])
  
  # Create boundary polygon
  boundary <- sf::st_as_sfc(site_poly$boundary_2157, crs = 2157)
  
  list(
    site_name = site,
    full_name = site_mapping[[site]],
    tower_x = coords[1],
    tower_y = coords[2],
    boundary_poly = boundary,
    raw_polygon = site_poly
  )
}

#------------------------------#
#  PBLH Calculation & QC       #
#------------------------------#

# Function to process each dataset
process_dataset <- function(df) {
  # Calculate PBLH and measurement height
  
  processed_df <- df %>%
    rowwise() %>%
    mutate(
      pblh = boundary_layer_height(
        MO_LENGTH_fluxnet, 
        u_star, 
        `w/ts_cov`, 
        ifelse(is.na(latitude_meta), 52.0, latitude_meta), 
        `(z-d)/L`, 
        air_temperature
      ),
      zm = master_sonic_height_meta - displacement_height_meta
    ) %>%
    ungroup() %>%
    filter(!is.na(zm), !is.na(MO_LENGTH_fluxnet), !is.na(roughness_length_meta), !is.na(pblh), !is.na(V_SIGMA_fluxnet), !is.na(u_star), !is.na(wind_dir))
  
  
  return(processed_df)
}

#------------------------------------------#
#  Define compute_footprint_ratio function #
#------------------------------------------#
compute_footprint_ratio <- function(footprint_result, boundary_poly) {
  # Ensure that the footprint result has the expected elements
  if (is.null(footprint_result$f_2d) ||
      is.null(footprint_result$x_2d) ||
      is.null(footprint_result$y_2d)) {
    return(0.0)
  }
  
  # Convert x, y, and flux matrices to vectors
  x_vec <- as.vector(footprint_result$x_2d)
  y_vec <- as.vector(footprint_result$y_2d)
  flux_vec <- as.vector(footprint_result$f_2d)
  
  # Create an sf object of points from the x,y coordinates (CRS = EPSG:2157)
  pts <- sf::st_as_sf(data.frame(x = x_vec, y = y_vec), coords = c("x", "y"), crs = 2157)
  
  # Determine which points are inside the boundary polygon
  inside <- st_within(pts, boundary_poly, sparse = FALSE)[, 1]
  
  # Calculate the total flux and the flux inside the polygon.
  total_flux <- sum(flux_vec, na.rm = TRUE)
  inside_flux <- sum(flux_vec[inside], na.rm = TRUE)
  
  # Compute and return the ratio (if total flux is > 0)
  ratio <- ifelse(total_flux > 0, inside_flux / total_flux, 0.0)
  return(ratio)
}

# Use calc_footprint_FFP_mod, translated from the Kljun FFP Python
add_footprint_ratio <- function(data, boundary_poly, tower_x, tower_y) {
  results <- data %>%
    rowwise() %>%
    mutate(
      footprint = list(
        calc_footprint_FFP_mod(
          zm = zm,
          z0 = roughness_length_meta,
          umean = wind_speed,
          h = pblh,
          ol = MO_LENGTH_fluxnet,
          sigmav = V_SIGMA_fluxnet,
          ustar = u_star,
          wind_dir = wind_dir,
          tower_x = tower_x,
          tower_y = tower_y
        )
      ),
      fp_ratio = compute_footprint_ratio(footprint, boundary_poly)
    ) %>%
    ungroup()
  
  data <- data %>%
    mutate(
      footprint_ratio = results$fp_ratio,
      Flux_QC_Footprint_Grid_50 = as.integer(footprint_ratio >= 0.5),
      Flux_QC_Footprint_Grid_60 = as.integer(footprint_ratio >= 0.6),
      Flux_QC_Footprint_Grid_70 = as.integer(footprint_ratio >= 0.7),
      Flux_QC_Footprint_Grid_80 = as.integer(footprint_ratio >= 0.8),
      Flux_QC_Footprint_Grid_90 = as.integer(footprint_ratio >= 0.9)
    )
  
  return(list(data = data, results = results))
}
