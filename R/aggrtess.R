#' Sum of point pattern marks over tessellation
#'
#' @param X point pattern of class `ppp`.
#' @param tess tessellation of class `tess`.
#'
#' @return Matrix of counts with one column for each tile in tessellation and
#' one row for each mark variable.
#' @export
#'
#' @examples
#' set.seed(42)
#' frb_centers <- rthin(frb, 0.01)
#' frb_tess <- augmented_dirichlet(frb_centers)
#' counts <- markquadratcount(frb, frb_tess)
#' marks(frb_tess) <- t(counts)
#' if(interactive()){
#'   plot(frb_tess, do.col = TRUE)
#' }
markquadratcount <- function(X, tess){
  Xsplit <- split(X, f = tess)
  counts <- sapply(Xsplit, function(x) colSums(as.data.frame(marks(x))))
  counts <- matrix(counts, ncol = length(Xsplit))
  return(counts)
}

validate_points_obj <- function(X){
  # Ensure X is a marked ppp
  stopifnot(is.ppp(X))
  if(!is.marked(X)){
    X$marks <- data.frame(count = rep(1, npoints(X)))
    X$markformat <- "dataframe"
  }
  return(X)
}

validate_target <- function(X, target){
  m <- marks(X, drop = FALSE)
  stopifnot(is.numeric(target) && ncol(m)==length(target))
  return(target)
}

#' Initialize aggregation tessellation
#'
#' @param X Point pattern with n_p points and n_m columns of marks
#' @param centers Point pattern of n_c tile centers or a string (see details).
#' @param target vector of length n_m
#' @param ... ignored.
#' @param adjust adjustment of target. One or n_m numbers between 0 and 1.
#'
#' @return An aggregation tessellation of class `aggr_tess` which is a list with components:
#'   - X: Original marked point pattern
#'   - centers: Point pattern of n_c tile centers
#'   - cells: tessellation with n_c tiles augmented with neighbour info
#'   - counts: n_m by n_c matrix with aggregates of mark values in each tile
#' Note: The aggregation tessellation does **not** necessarily satisfy the aggregation
#' target at this point.
#' @export
#'
#' @examples
#' aggrtess(frb, centers = "rthin", target = c(50, 100), )
aggrtess <- function(X, centers, target, ..., adjust = 1){
  X <- validate_points_obj(X)
  target <- validate_target(X, target)

  # Interpret centers
  if(!inherits(centers, "ppp")){
    stopifnot(is.character(centers))
    centers <- make_aggrtess_centers(X, centers, target, ..., adjust = adjust)
  }

  tessel <- augmented_dirichlet(centers)
  counts <- markquadratcount(X, tess = tessel)

  out <- list(X = X, tessel = tessel, counts = counts, target = target)
  class(out) <- "aggrtess"
  return(out)
}

#' Title
#'
#' @param X Point pattern with n_p points and n_m columns of marks
#' @param centers String (see details).
#' @param target vector of length n_m
#' @param ... ignored.
#' @param adjust adjustment of target. One or n_m numbers between 0 and 1.
#'
#' @return A point pattern of class `ppp`.
#' @export
#'
#' @examples
#' make_aggrtess_centers(frb, "rthin", c(50, 100))
make_aggrtess_centers <- function(X, centers, target, ..., adjust = 1){
  # Validate
  m <- marks(X, drop = FALSE)
  stopifnot(length(adjust)==1 | length(adjust)==length(target))
  # Decide on number of tiles
  ntiles <- min(colSums(m)/(adjust*target))
  nX <- npoints(X)
  ntiles <- min(ntiles, nX)
  # Do easy cases
  if(centers == "rthin"){
    return(rthin(X, P = ntiles/nX))
  }
  if(centers == "index_rand"){
    return(X[sample(nX, ntiles)])
  }
  if(centers == "index_shift"){
    step <- floor(nX/ntiles)
    start <- sample(step,1)
    return(X[seq(start, nX, by = step)])
  }

  ### Now starts grid-like centers ###
  W <- Frame(X)
  # Grow W to square
  dx <- diff(W$xrange)
  dy <- diff(W$yrange)
  if(dx<dy){
    W_ext <- grow.rectangle(W, xmargin = (dy-dx)/2, ymargin = 0)
  } else{
    W_ext <- grow.rectangle(W, xmargin = 0, ymargin = (dx-dy)/2)
  }
  ntiles_ext <- ntiles*area(W_ext)/area(W)
  nx_ext <- ceiling(sqrt(ntiles_ext))
  if(centers == "rectgrid"){
    cent_ext <- as.ppp(gridcenters(W_ext, nx = nx_ext, ny = nx_ext), W = W_ext)
    return(cent_ext[Window(X)])
  }
  if(centers == "hexgrid"){
    cent_ext <- hexgrid(W_ext, s = sqrt(2/(3*sqrt(3)))*sidelengths(W_ext)[1]/nx_ext)
    return(cent_ext[Window(X)])
  }

  #### From this point remains only quantess ####
  stop("Not implemented yet!")
}
