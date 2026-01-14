#' Create a bounding box
#'
#' This function creates a bounding box from either cardinal coordinates
#' (north, south, east and west) or center coordinates and dimensions
#' (center longitude and latitude, and size of width and height).
#'
#' @param north The northernmost latitude.
#' @param south The southernmost latitude.
#' @param east The easternmost longitude.
#' @param west The westernmost longitude.
#' @param center_lon The center longitude.
#' @param center_lat The center latitude.
#' @param width The width of the bounding box.
#' @param height The height of the bounding box.
#' @return A 2x2 matrix representing the bounding box.
#' @export
#' @examples
#' # Create a bounding box from cardinal coordinates
#' create_bbox(
#'   north =  35.1899,
#'   south =  35.1399,
#'   east  = 136.9524,
#'   west  = 136.8524
#' )
#'
#' # Create a bounding box from center coordinates and dimensions
#' create_bbox(
#'   center_lon = 136.9024,
#'   center_lat =  35.1649,
#'   width      = 0.10,
#'   height     = 0.05
#' )
create_bbox <- function(north = NULL, south = NULL,
                        east = NULL, west = NULL,
                        center_lon = NULL, center_lat = NULL,
                        width = NULL, height = NULL) {
  if (!is.null(north) && !is.null(south) && !is.null(east) && !is.null(west)) {
    validate_cardinal_points(north, south, east, west)
  } else if (!is.null(center_lon) && !is.null(center_lat) &&
             !is.null(width) && !is.null(height)) {
    validate_dimensions(width, height)
    west <- center_lon - width / 2
    east <- center_lon + width / 2
    south <- center_lat - height / 2
    north <- center_lat + height / 2
  } else {
    stop("invalid argument provided")
  }

  bbox <- matrix(c(west, south, east, north), nrow = 2)
  dimnames(bbox) <- list(c("x", "y"), c("min", "max"))
  return(bbox)
}

validate_cardinal_points <- function(north, south, east, west) {
  if (north <= south) stop("`north` must be greater than `south`")
  if (east <= west) stop("`east` must be greater than `west`")
}

validate_dimensions <- function(width, height) {
  if (width <= 0 || height <= 0) stop("`width` and `height` must be positive")
}



#' Create a bounding box from center coordinates and distance in km
#'
#' @description
#' Calculates a rectangular bounding box based on a central coordinate and
#' physical dimensions in kilometers. Uses JGD2011 Japan Plane Rectangular
#' Coordinate System for accurate distance calculation in Japan, and returns
#' the bounding box in WGS84 as a 2x2 matrix.
#'
#' @param center_lon Center longitude (degrees)
#' @param center_lat Center latitude (degrees)
#' @param width_km Width of the rectangle in km (East-West)
#' @param height_km Height of the rectangle in km (North-South)
#' @param jgd2011_zone JGD2011 Plane Rectangular CS zone number (1-19).
#'   Default is 7 (covers Aichi, Mie, Gifu, etc.). See details for zone coverage.
#' @return A 2x2 matrix representing the bounding box
#' @details
#' JGD2011 Plane Rectangular Coordinate System zones for Japan:
#' \itemize{
#'   \item Zone 1 (EPSG:6669): Nagasaki, Kagoshima (islands)
#'   \item Zone 2 (EPSG:6670): Fukuoka, Saga, Nagasaki, Kumamoto, Oita, Miyazaki, Kagoshima
#'   \item Zone 3 (EPSG:6671): Yamaguchi, Shimane, Hiroshima
#'   \item Zone 4 (EPSG:6672): Okayama, Tottori, Hyogo (part), Kagawa, Tokushima, Ehime, Kochi
#'   \item Zone 5 (EPSG:6673): Hyogo, Osaka, Kyoto, Wakayama, Nara, Shiga
#'   \item Zone 6 (EPSG:6674): Shizuoka, Yamanashi, Nagano, Gifu (part)
#'   \item Zone 7 (EPSG:6675): Aichi, Mie, Gifu, Ishikawa, Fukui, Toyama
#'   \item Zone 8 (EPSG:6676): Niigata, Nagano (part), Gunma, Tochigi
#'   \item Zone 9 (EPSG:6677): Tokyo, Chiba, Saitama, Kanagawa, Ibaraki
#'   \item Zone 10-19: Northern regions and islands
#' }
#' @examples
#' # Create a 6.5 km x 4.5 km bounding box centered around Nagoya (Zone 7)
#' bbox <- create_bbox_km_cartesian(
#'   center_lon = 136.9,
#'   center_lat = 35.17,
#'   width_km   = 6.5,
#'   height_km  = 4.5
#' )
#' @export
create_bbox_km_cartesian <- function(center_lon, center_lat, width_km, height_km,
                                     jgd2011_zone = 7) {

  if (any(sapply(list(center_lon, center_lat, width_km, height_km), is.null))) {
    stop("All arguments must be provided")
  }
  if (width_km <= 0 || height_km <= 0) stop("width_km and height_km must be positive")
  if (jgd2011_zone < 1 || jgd2011_zone > 19) stop("jgd2011_zone must be between 1 and 19")

  # JGD2011 Plane Rectangular CS: EPSG:6669 (zone 1) to EPSG:6687 (zone 19)
  jgd2011_epsg <- 6668 + jgd2011_zone

  # Transform center point to JGD2011 Cartesian coordinates
  center_wgs84 <- sf::st_sfc(sf::st_point(c(center_lon, center_lat)), crs = 4326)
  center_jgd <- sf::st_transform(center_wgs84, crs = jgd2011_epsg)

  coords <- sf::st_coordinates(center_jgd)
  x0 <- coords[1]
  y0 <- coords[2]

  # Calculate offsets in meters
  dx <- width_km * 1000 / 2
  dy <- height_km * 1000 / 2

  # Calculate corner coordinates in JGD2011
  west_jgd  <- x0 - dx
  east_jgd  <- x0 + dx
  south_jgd <- y0 - dy
  north_jgd <- y0 + dy

  # Transform corners back to WGS84
  corners_jgd <- sf::st_sfc(
    sf::st_point(c(west_jgd, south_jgd)),
    sf::st_point(c(east_jgd, north_jgd)),
    crs = jgd2011_epsg
  )
  corners_wgs84 <- sf::st_transform(corners_jgd, 4326)
  corner_coords <- sf::st_coordinates(corners_wgs84)

  west  <- corner_coords[1, 1]
  south <- corner_coords[1, 2]
  east  <- corner_coords[2, 1]
  north <- corner_coords[2, 2]

  bbox <- matrix(c(west, south, east, north), nrow = 2)
  dimnames(bbox) <- list(c("x", "y"), c("min", "max"))

  return(bbox)
}
