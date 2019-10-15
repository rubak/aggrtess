#' Script to download polygon outline and all addresses in the Danish
#' 'kommune' Frederiksberg and attach artificial information about the
#' number of households and the number of individuals associated with
#' each address.

library(here)
library(geojsonsf)
library(sf)
library(maptools)
library(spatstat)

# Download 'adgangsadresser' in Frederiksberg and keep x,y coords.
www_pts <- "dawa.aws.dk/adgangsadresser?kommunekode=0147&format=csv&srid=25832&struktur=mini"
tmpfile <- tempfile(fileext = "csv")
download.file(www_pts, destfile = tmpfile)
pts <- read.csv(tmpfile)
unlink(tmpfile)
pts <- pts[, c("x", "y")]

# Get polygon data
www_poly <- "http://dawa.aws.dk/kommuner/0147?format=geojson&srid=25832"
poly <- geojsonsf::geojson_sfc(www_poly)

# Convert to statstat format
W <- as(sf::as_Spatial(poly), "owin")
frb <- as.ppp(pts, W = W)

## Generate artificial household and individuals data:
set.seed(42)
## First one household at all addresses
households <- rep(1, npoints(frb))
## Sample some households and assign them a larger number of households
i <- sample(npoints(frb), npoints(frb)/10)
households[i] <- 2 + rpois(length(i), 1)
## Generate an average number of people for each household and multiply and round
individuals <- runif(npoints(frb), 1, 10)
individuals[individuals>4] <- 1
individuals <- round(households * individuals)

# Attach data as marks to points
marks(frb) <- data.frame(households, individuals)

save(frb, file = here::here("data", "frb.rda"),
     version = 2, compress = "bzip2")
