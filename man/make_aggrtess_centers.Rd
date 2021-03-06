% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggrtess.R
\name{make_aggrtess_centers}
\alias{make_aggrtess_centers}
\title{Generate tile centers for an aggregation tessellation}
\usage{
make_aggrtess_centers(
  X,
  centers = c("rectgridsample", "hexgridsample", "index_rand", "index_shift", "rthin",
    "rectgridcenter", "hexgridcenter"),
  target,
  ...,
  adjust = 1
)
}
\arguments{
\item{X}{Point pattern with n_p points and n_m columns of marks}

\item{centers}{String to specify the type of center distribution (see details).}

\item{target}{vector of length n_m}

\item{...}{ignored.}

\item{adjust}{adjustment of target. One or n_m positive numbers. E.g.
\code{adjust = 2} doubles the target count in each tile, so the number of tiles
(centers) is approximately halved.}
}
\value{
A point pattern of class \code{ppp}.
}
\description{
Generate tile centers for an aggregation tessellation
}
\details{
The arguments \code{target} and \code{adjust} are used to determine
approximately how many tile centers to generate. Argument \code{centers} specifies
how to generate the centers. Current valid choices are:
\code{"rectgridsample"}, \code{"hexgridsample"}, \code{"index_rand"}, \code{"index_shift"},
\code{"rthin"}, \code{"rectgridcenter"}, \code{"hexgridcenter"}.
Methods starting with \code{"rectgrid"} and \code{"hexgrid"} start by making a rectangular
or hexagonal grid over the window of \code{X} and the either use the grid centers
or sample a point in each grid cell. Method \code{"rthin"} thins \code{X} and uses the
thinned pattern as tile centers (number of centers is then random and
subsequent identical calls to \code{make_aggrtess_centers} typically have a
different number of centers). Method \code{"index_rand"} samples a subset of \code{X}
of fixed size and \code{"index_shift"} shifts the starting index randomly and then
picks every n'th point of \code{X} from that starting index.
}
\examples{
x <- make_aggrtess_centers(frb, "hexgridcenter", target = c(50, 100))
x
if(interactive()){
plot(x)
}
}
