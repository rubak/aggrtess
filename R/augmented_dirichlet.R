#' Dirichlet tessellation augmented with neighbour info
#'
#' @param X point pattern
#'
#' @return Dirichlet tessellation of class `tess` with additional field `neigh`
#' @export
#' @importFrom deldir tile.list
#' @import spatstat
#'
#' @examples
#' set.seed(42)
#' frb_small <- rthin(frb, 0.01)
#' d <- augmented_dirichlet(frb_small)
#' head(d$neigh)
augmented_dirichlet <- function(X) {
  stopifnot(is.ppp(X))
  X <- unique(X, rule="deldir", warn=TRUE)
  nX <- npoints(X)
  w <- X$window
  if(nX == 0) return(NULL)
  if(nX == 1) return(as.tess(w))
  dd <- safedeldir(X)
  if(is.null(dd)) return(NULL)
  pp <- lapply(deldir::tile.list(dd), function(z) { owin(poly=z[c("x","y")]) })
  if(length(pp) == npoints(X))
    names(pp) <- seq_len(npoints(X))
  dir <- tess(tiles=pp, window=as.rectangle(w))
  if(w$type != "rectangle")
    dir <- intersect.tess(dir, w)
  # dir$dirsgs <- dd$dirsgs
  # Epsilon for numerical comparison below
  eps <- sqrt(.Machine$double.eps)
  # Making the neighbour list
  neigh <- vector("list", length = dir$n)
  tilelist <- tiles(dir)
  vertices_list <- lapply(tilelist, vertices.owin)
  vertices_list <- lapply(vertices_list, function(v){ppp(v$x, v$y, window = w, check = FALSE)})
  for(i in seq_len(dir$n)){
    ## All neighbours from deldir (infinite border tiles)
    all_neigh <- c(dd$delsgs$ind1[dd$delsgs$ind2==i],
                   dd$delsgs$ind2[dd$delsgs$ind1==i])
    ## The remainder keeps only neighbour tiles that share a vertex with tile i:
    true_neigh <- sapply(vertices_list[all_neigh], function(xx){min(nncross.ppp(vertices_list[[i]], xx))}) < eps
    neigh[[i]] <- all_neigh[true_neigh]
  }
  dir$neigh <- neigh
  return(dir)
}
