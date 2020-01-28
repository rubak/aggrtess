#' Conversion of aggrtess to sf format
#'
#' @param x tessellation
#' @param crs crs to use with sf
#'
#' @return Simple features object with the polygons and their counts.
#' @export
#' @importFrom methods as
#'
#' @examples
#'
#' set.seed(42)
#' x <- aggrtess(frb, centers = "index_shift", target = c(50, 100))
#' x_sf <- aggrtess2sf(x)
#' if(interactive()){
#'   mapview::mapview(x_sf, zcol = c("households", "individuals"))
#' }
aggrtess2sf <- function(x, crs = 25832){
  if(! requireNamespace("maptools") || ! requireNamespace("sf")){
    stop("Packages maptools and sf are needed. Please install.")
  }
  poly <- lapply(tiles(x$tessel), function(y) sf::st_as_sf(as(y, "SpatialPolygons"), crs = crs))
  poly <- do.call(rbind, poly)
  sf::st_crs(poly) <- crs
  counts <- t(x$counts)
  colnames(counts) <- colnames(marks(x$X))
  dat <- data.frame(id = 1:nrow(counts), counts, poly)
  return(sf::st_sf(dat))
}
