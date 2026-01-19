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
#' @param center_lon The center longitude (WGS84).
#' @param center_lat The center latitude (WGS84).
#' @param width The width of the bounding box.
#' @param height The height of the bounding box.
#' @param unit Unit for width and height: \code{"degree"} (default) or
#'   \code{"km"}. When \code{"km"} is specified, width and height are
#'   interpreted as kilometers and a projected CRS is used for calculation.
#' @param crs EPSG code of a projected coordinate system to use when
#'   \code{unit = "km"}. Default is 6675 (JGD2011 Zone 7). Ignored when
#'   \code{unit = "degree"}.
#' @return A 2x2 matrix representing the bounding box in WGS84 coordinates.
#' @details
#' When \code{unit = "km"}, common projected CRS options for Japan are:
#' \itemize{
#'   \item EPSG:6669 (Zone 1): Nagasaki, Kagoshima (islands)
#'   \item EPSG:6670 (Zone 2): Fukuoka, Saga, Kumamoto, Oita, Miyazaki, Kagoshima
#'   \item EPSG:6671 (Zone 3): Yamaguchi, Shimane, Hiroshima
#'   \item EPSG:6672 (Zone 4): Okayama, Tottori, Kagawa, Tokushima, Ehime, Kochi
#'   \item EPSG:6673 (Zone 5): Hyogo, Osaka, Kyoto, Wakayama, Nara, Shiga
#'   \item EPSG:6674 (Zone 6): Shizuoka, Yamanashi, Nagano, Gifu (part)
#'   \item EPSG:6675 (Zone 7): Aichi, Mie, Gifu, Ishikawa, Fukui, Toyama
#'   \item EPSG:6676 (Zone 8): Niigata, Nagano (part), Gunma, Tochigi
#'   \item EPSG:6677 (Zone 9): Tokyo, Chiba, Saitama, Kanagawa, Ibaraki
#' }
#' UTM zones can also be used (e.g., EPSG:32654 for UTM Zone 54N).
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
#' # Create a bounding box from center coordinates and dimensions (degrees)
#' create_bbox(
#'   center_lon = 136.9024,
#'   center_lat =  35.1649,
#'   width      = 0.10,
#'   height     = 0.05
#' )
#'
#' # Create a bounding box with dimensions in kilometers
#' create_bbox(
#'   center_lon = 136.9,
#'   center_lat = 35.17,
#'   width      = 6.5,
#'   height     = 4.5,
#'   unit       = "km"
#' )
#'
#' # Using a different CRS (Tokyo area with JGD2011 Zone 9)
#' create_bbox(
#'   center_lon = 139.7,
#'   center_lat = 35.68,
#'   width      = 5.0,
#'   height     = 5.0,
#'   unit       = "km",
#'   crs        = 6677
#' )
create_bbox <- function(north = NULL, south = NULL,
                        east = NULL, west = NULL,
                        center_lon = NULL, center_lat = NULL,
                        width = NULL, height = NULL,
                        unit = c("degree", "km"),
                        crs = 6675) {
  unit <- match.arg(unit)

  if (!is.null(north) && !is.null(south) && !is.null(east) && !is.null(west)) {
    validate_cardinal_points(north, south, east, west)
  } else if (!is.null(center_lon) && !is.null(center_lat) &&
             !is.null(width) && !is.null(height)) {
    validate_dimensions(width, height)

    if (unit == "degree") {
      west <- center_lon - width / 2
      east <- center_lon + width / 2
      south <- center_lat - height / 2
      north <- center_lat + height / 2
    } else {
      # unit == "km"
      coords <- calculate_bbox_from_km(center_lon, center_lat, width, height, crs)
      west <- coords$west
      east <- coords$east
      south <- coords$south
      north <- coords$north
    }
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

calculate_bbox_from_km <- function(center_lon, center_lat, width_km, height_km, crs) {
  # Transform center point to projected coordinates
  center_wgs84 <- sf::st_sfc(sf::st_point(c(center_lon, center_lat)), crs = 4326)
  center_projected <- sf::st_transform(center_wgs84, crs = crs)

  coords <- sf::st_coordinates(center_projected)
  x0 <- coords[1]
  y0 <- coords[2]

  # Calculate offsets in meters
  dx <- width_km * 1000 / 2
  dy <- height_km * 1000 / 2

  # Calculate corner coordinates in projected CRS
  west_proj  <- x0 - dx
  east_proj  <- x0 + dx
  south_proj <- y0 - dy
  north_proj <- y0 + dy

  # Transform corners back to WGS84
  corners_proj <- sf::st_sfc(
    sf::st_point(c(west_proj, south_proj)),
    sf::st_point(c(east_proj, north_proj)),
    crs = crs
  )
  corners_wgs84 <- sf::st_transform(corners_proj, 4326)
  corner_coords <- sf::st_coordinates(corners_wgs84)

  list(
    west  = corner_coords[1, 1],
    south = corner_coords[1, 2],
    east  = corner_coords[2, 1],
    north = corner_coords[2, 2]
  )
}
