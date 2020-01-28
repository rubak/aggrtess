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
#' @param ... passed to `make_aggrtess_centers` when `target` is a string..
#'
#' @return An aggregation tessellation of class `aggr_tess` which is a list with
#' the main components being:
#'   - X: Original marked point pattern
#'   - centers: Point pattern of n_c tile centers
#'   - tessel: tessellation with n_c tiles augmented with neighbour info
#'   - counts: n_m by n_c matrix with aggregates of mark values in each tile
#' Note: The aggregation tessellation does **not** necessarily satisfy the aggregation
#' target at this point.
#' @export
#'
#' @examples
#' x <- aggrtess(frb, centers = "rectgridsample", target = c(50, 100))
#' x
aggrtess <- function(X, centers, target, ...){
  X <- validate_points_obj(X)
  target <- validate_target(X, target)

  # Interpret centers
  if(!inherits(centers, "ppp")){
    stopifnot(is.character(centers))
    centers <- make_aggrtess_centers(X, centers, target, ...)
  }
  # Make sure centers are unmarked
  centers <- unmark(centers)

  tessel <- augmented_dirichlet(centers)
  counts <- markquadratcount(X, tess = tessel)
  deviation <- (counts - target)
  clust <- calc_cluster_size(counts, tessel$neigh)
  excess <- clust$total - target * matrix(clust$n, nrow(clust$total), ncol(clust$total), byrow = TRUE)

  out <- list(X = X, centers = centers, tessel = tessel, target = target,
              counts = counts, deviation = deviation, excess = excess,
              clustersize = clust)
  class(out) <- "aggrtess"
  return(out)
}

#' @export
print.aggrtess <- function(x, ...){
  cat("Aggregation tessellation with:",
      paste(npoints(x$X), "points split into", x$tessel$n, "tiles."),
      paste("Marks being aggregated:", paste(names(marks(x$X)), collapse = ",")),
      paste("Aggregation target:", paste(x$target, collapse = ",")),
      sep = "\n")
}

#' Plot aggregation tessellation of class aggrtess
#'
#' @param x aggregation tessellation object of class `aggrtess`
#' @param ... passed to `plot.tess` in spatstat.
#' @param what string to determine what to plot. Either `"counts"`,
#' `"deviation"` or `"excess"`.
#' @param collapse function to use to collapse multiple values to a single value
#' if desired. Typically `sum` or `mean`.
#' @param relative logical to use relative values. Divides by number of tiles in
#' cluster for `excess` and divides by `x$target` otherwise.
#' @param do.col color each tile according to associated value. Defaults to `TRUE`.
#' @param do.labels annotate plot with tile labels (ids). Defaults to `TRUE`.
#' @param main main title of plot. Defaults to empty string.
#'
#' @return the marked tessellation of class `tess` which is plotted (invisibly).
#' @importFrom spatstat.utils resolve.defaults
#' @export
#'
#' @examples
#' if(interactive()){
#' x <- aggrtess(frb, centers = "rectgridsample", target = c(50, 100))
#' plot(x)
#' }
plot.aggrtess <- function(x, ..., what = "counts", collapse = NULL,
                          relative = FALSE, do.col = TRUE, do.labels = TRUE, main = ""){
  match.arg(what, c("counts", "deviation", "excess"))
  denom <- 1 # Divisor for relative = TRUE
  if(relative){
    denom <- x$target
    if(what == "excess"){
      size <- x$clustersize
      denom <- matrix(size$n, nrow(size$total), ncol(size$total), byrow = TRUE)
    }
  }
  m <- x[[what]]/denom
  if(!is.null(collapse)){
    stopifnot(is.function(collapse))
    m <- apply(m, 2, collapse)
  }
  m <- as.data.frame(t(m))
  names(m) <- if(is.null(collapse)){
    names(marks(x$X))
  } else{
      as.character(substitute(collapse))
  }
  tmp <- x$tessel %mark% m
  dots <- spatstat.utils::resolve.defaults(
    ..., list(do.col = do.col, do.labels = do.labels, main = main))
  do.call("plot.tess", append(list(x=tmp), dots))
  return(tmp)
}

#' Find deviation from target for aggegation tessellation
#'
#' @param x object of class `aggrtess`
#' @param ... ignored
#' @param which character describing the desired tiles: `"all"` (all `n` tiles),
#' `"bad"` (0-`n` tiles not satisfying target), or `"worst"` (single tile).
#' @param relative logical to use relative deviations from target.
#'
#' @return `data.frame` with columns `id`, `min` and one additional column for
#' each mark. Number of rows is between 0 and `n` depending on `which`.
#' @export
#'
#' @examples
#' x <- aggrtess(frb, centers = "rectgridsample", target = c(50, 100))
#' target_deviation(x)
target_deviation <- function(x, ..., which = c("all", "bad", "worst"), relative = TRUE){
  which <- match.arg(which)
  dev <- x$deviation
  target <- x$target
  if(relative){
    dev <- dev/target
  }
  out <- as.data.frame(t(dev))
  names(out) <- names(marks(x$X))
  out_min <- apply(out, 1, min)
  # out_mean <- apply(out, 1, mean)
  # out_max <- apply(out, 1, max)
  out <- cbind(id = seq_len(nrow(out)), out, min = out_min)
  if(which == "bad"){
    out <- out[which(out$min<0),]
  }
  if(which == "worst"){
    out <- out[which.min(out$min),]
  }
  return(out)
}

calc_cluster_size <- function(counts, neigh){
  size <- list(n=rep(0, ncol(counts)), total = 0*counts)
  for(i in seq_along(neigh)){
    size$n[i] <- length(neigh[[i]])+1
    size$total[,i] <- rowSums(counts[,c(i, neigh[[i]]),drop=FALSE])
  }
  return(size)
}


#' Generate tile centers for an aggregation tessellation
#'
#' @param X Point pattern with n_p points and n_m columns of marks
#' @param centers String to specify the type of center distribution (see details).
#' @param target vector of length n_m
#' @param ... ignored.
#' @param adjust adjustment of target. One or n_m positive numbers. E.g.
#' `adjust = 2` doubles the target count in each tile, so the number of tiles
#' (centers) is approximately halved.
#'
#' @return A point pattern of class `ppp`.
#' @details The arguments `target` and `adjust` are used to determine
#' approximately how many tile centers to generate. Argument `centers` specifies
#' how to generate the centers. Current valid choices are:
#' `"rectgridsample"`, `"hexgridsample"`, `"index_rand"`, `"index_shift"`,
#' `"rthin"`, `"rectgridcenter"`, `"hexgridcenter"`.
#' Methods starting with `"rectgrid"` and `"hexgrid"` start by making a rectangular
#' or hexagonal grid over the window of `X` and the either use the grid centers
#' or sample a point in each grid cell. Method `"rthin"` thins `X` and uses the
#' thinned pattern as tile centers (number of centers is then random and
#' subsequent identical calls to `make_aggrtess_centers` typically have a
#' different number of centers). Method `"index_rand"` samples a subset of `X`
#' of fixed size and `"index_shift"` shifts the starting index randomly and then
#' picks every n'th point of `X` from that starting index.
#' @export
#'
#' @examples
#' x <- make_aggrtess_centers(frb, "hexgridcenter", target = c(50, 100))
#' x
#' if(interactive()){
#' plot(x)
#' }
make_aggrtess_centers <-
  function(X,
           centers = c("rectgridsample", "hexgridsample", "index_rand", "index_shift", "rthin", "rectgridcenter", "hexgridcenter"),
           target, ..., adjust = 1){
  # Validate
  centers <- match.arg(centers)
  m <- marks(X, drop = FALSE)
  X <- unmark(X)
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
  if(centers == "rectgridcenter"){
    cent_ext <- as.ppp(gridcenters(W_ext, nx = nx_ext, ny = nx_ext), W = W_ext)
    return(cent_ext[Window(X)])
  }
  if(centers == "hexgridcenter"){
    cent_ext <- hexgrid(W_ext, s = sqrt(2/(3*sqrt(3)))*sidelengths(W_ext)[1]/nx_ext)
    return(cent_ext[Window(X)])
  }
  if(centers %in% c("rectgridsample", "hexgridsample")){
    if(centers == "rectgridsample"){
      grid_ext <- quadrats(W_ext, nx = nx_ext, ny = nx_ext)
    } else{
      grid_ext <- hextess(W_ext, s = sqrt(2/(3*sqrt(3)))*sidelengths(W_ext)[1]/nx_ext)
    }
    X_split <- split(X, grid_ext, drop = TRUE)
    names(X_split) <- NULL
    X_sample <- solapply(X_split, function(x) x[sample(npoints(x), 1)])
    cent <- superimpose(X_sample, W = Window(X), check = FALSE)
    return(cent)
  }

  #### From this point remains only quantess ####
  stop("Not implemented yet!")
}

#' Delete tile in aggregation tessellation
#'
#' @param x Aggregation tessellation of class `aggrtess`.
#' @param id Integer indicating which tile to delete
#'
#' @return Aggregation tessellation. Either with one less tile than the input
#' (if split was successful) or with the same number of tiles (if split didn't succeed in `max_tries` attempts).
#' @export
#'
#' @examples
#' x <- aggrtess(frb, centers = "hexgridcenter", target = c(50, 100))
#' y <- delete_tile(x, 4)
delete_tile <- function(x, id){
  stopifnot(is.numeric(id) && all(id>=1) && all(id<=npoints(x$centers)))
  return(aggrtess(x$X, x$centers[-id], x$target))
}

#' Prune aggregation tessellation to remove all tiles below the target
#'
#' @param x Aggregation tessellation of class `aggrtess`.
#' @param max_tiles Maximal number of tiles to delete (to ensure the computation ends in a resonable time if started in error on a very bad/big example.)
#' @param ... Ignored.
#' @param verbose Logical to print diagnostics.
#'
#' @return Aggregation tessellation.
#' @export
#'
#' @examples
#' X <- frb[seq(1, npoints(frb), by = 10)]
#' x <- aggrtess(X, centers = "hexgridcenter", target = c(50, 100))
#' y <- prune_aggrtess(x)
prune_aggrtess <- function(x, max_tiles = 100, ..., verbose = FALSE){
  iter <- 0
  worst <- target_deviation(x, which = "worst")
  while(worst$min<0 & iter<max_tiles){
    iter <- iter+1
    if(verbose){
      nbad <- nrow(target_deviation(x, which = "bad"))
      print("Iter.", iter, "Currently", nbad, "tiles are bad.")
    }
    x <- delete_tile(x, id = worst$id)
    worst <- target_deviation(x, which = "worst")
  }
  return(x)
}

#' Split tile in aggregation tessellation
#'
#' @param x Aggregation tessellation of class `aggrtess`.
#' @param id Integer id of the tile to try to split
#' @param ... Ignored
#' @param max_tries Maximal number of times to try to split tile before giving up.
#' @param verbose Logical to print diagnostics.
#'
#' @return Aggregation tessellation. Either with one more tile than the input
#' (if split was successful) or with the same number of tiles (if split didn't succeed in `max_tries` attempts).
#' @export
#'
#' @examples
#' X <- frb[seq(1, npoints(frb), by = 10)]
#' x <- aggrtess(X, centers = "hexgridcenter", target = c(50, 100))
#' x <- prune_aggrtess(x)
#' set.seed(42) # Reproducibility
#' y <- split_tile(x, 6)
split_tile <- function(x, id, ..., max_tries = 100, verbose = FALSE){
  box <- Frame(x$X)
  Wlist <- tiles(x$tessel)
  # Window of tile `id` and the points in this tile
  W_id <- Wlist[[id]]
  X_id <- x$X[W_id]
  # Neighbours of tile `id`
  neigh <- x$tessel$neigh[[id]]
  # Cluster around tile `id`, i.e., tile `id` and its neighbours
  Wclust <- do.call(union.owin, Wlist[c(id,neigh)])
  Xclust <- x$X[Wclust]
  iter <- 0
  OK <- FALSE
  while(iter<max_tries & !OK){
    iter <- iter + 1
    # Propose two new tile centers inside tile `id`
    centers_id <- unmark(X_id[sample(npoints(X_id), 2)])
    # Combine proposed centers with neighbouring centers and check validity of cluster
    centersclust <- superimpose.ppp(x$centers[neigh], centers_id, W = Wclust, check = FALSE)
    tesselclust <- dirichlet(centersclust)
    countclust <- markquadratcount(Xclust, tess = tesselclust)
    OK <- all(countclust>x$target)
    # Need full check of validity below due to non-convexity at boundaries etc.
    if(OK){ # Potentially OK
      if(verbose){
        print(paste("Potential split of cell", id, "in iter.", iter, "/", max_tries))
      }
      x_prop <- aggrtess(x$X, superimpose(x$centers[-id], centers_id), x$target)
      OK <- target_deviation(x_prop, which = "worst")$min >= 0
      if(OK){
        if(verbose){
          print(paste("succesfully changed tile", id, "in iter.", iter, "/", max_tries))
        }
        x <- x_prop
      } else{
        if(verbose){
          print(paste("Sorry. No success for this one!"))
        }
      }
    }
  }
  attr(x, "changed") <- OK
  return(x)
}

#' Fill in tiles in aggregation tessellation by sequential splitting
#'
#' @param x Aggregation tessellation of class `aggrtess`.
#' @param ... Ignored
#' @param max_tiles Maximal number of tiles to try to split before giving up.
#' @param max_tries Maximal number of times to try to split each tile before giving up.
#' @param adjust adjustment of target. See details.
#' @param verbose Logical to print diagnostics.
#'
#' @details For each tile `i` the total aggregation count in the neighbourhood (tiles that share border with `i`) is calculated.
#' The size of the neighbourhood (including `i`) together with the target determines the minimum count required for the given neigbourhood size.
#' The difference between the actual count and the minimum required is called the excess. If the excess is larger than the target
#' the neighbourhood could in principle contain an extra tile and is eligible for a split. The one with the biggest excess is tried first.
#' When eligibility is considered an adjusted target (`adjust*target`) is considered rather than the raw target to avoid trying to split very difficult neighbourhoods
#' @return Aggregation tessellation.
#' @export
#'
#' @examples
#' X <- frb[seq(1, npoints(frb), by = 10)]
#' x <- aggrtess(X, centers = "hexgridcenter", target = c(50, 100))
#' x <- prune_aggrtess(x)
#' set.seed(42) # Reproducibility
#' y <- fill_aggrtess(x)
fill_aggrtess <- function(x, ..., max_tiles = 10, max_tries = 20, adjust = 1.2, verbose = FALSE){
  size <- x$clustersize
  excess <- size$total - x$target * matrix(size$n, nrow(size$total), ncol(size$total), byrow = TRUE)
  candidates <- which(apply(excess>adjust*x$target, 2, all))
  iter <- 0
  while(length(candidates)>0 & iter<max_tiles){
    excess_candidates <- excess[,candidates,drop=FALSE]
    nr <- nrow(excess_candidates)
    nc <- ncol(excess_candidates)
    avg_excess_candidates <- excess_candidates/matrix(size$n[candidates], nr, nc, byrow = TRUE)
    top <- col(avg_excess_candidates)[which.max(avg_excess_candidates)]
    iter <- iter+1
    ii <- candidates[top]
    if(verbose){
      print(paste("Iter.", iter, "/", max_tiles, "Trying tile", ii, "with", size$n[ii], "neighbours and excess", paste(excess[,ii], collapse = ",")))
    }
    tmp <- split_tile(x, ii, max_tries = max_tries, verbose = verbose)
    if(attr(tmp, "changed")){
      # print("Change detected!")
      x <- tmp
      size <- x$clustersize
      excess <- size$total - x$target * matrix(size$n, nrow(size$total), ncol(size$total), byrow = TRUE)
      candidates <- which(apply(excess>adjust*x$target, 2, all))
    } else{
      candidates <- candidates[-top]
    }
  }
  if(verbose && length(candidates)==0){
    print("No more candidate clusters for splitting.")
  }
  return(x)
}
