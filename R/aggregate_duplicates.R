utils::globalVariables(c(".SD", "id", ":="))

#' Aggregation of duplicated columns in a data.frame
#'
#' @param x data.frame with
#' @param idcols names of columns identifying the locations which should be checked for duplication. If `NULL` all other columns than `sumcols` are used. Can't be `NULL` if `sumcols` is.
#' @param sumcols names of columns to be aggregated for identical ids. If `NULL` all other columns than `idcols` are used. Can't be `NULL` if `idcols` is.
#'
#' @details **WARNING**: Factors are replaced by strings in the result.
#'
#' @return data.frame with columns as specified by `idcols` and `sumcols`.
#' Columns not (implicitly) given in these are dropped.
#' @export
#' @rawNamespace import(data.table, except = shift)
#' @importFrom utils globalVariables
#'
#' @examples
#' # Example with a single duplicated location
#' x <- data.frame(municipality = c("A", "B", "B"), x = c(1, 2, 2), y = c(1, 1, 1),
#'                 households = c(1, 1, 2), individuals = c(4, 3, 8))
#' # Three ways to achieve same aggregated version:
#' y1 <- aggregate_duplicates(x, idcols = c("municipality", "x", "y"),
#'                            sumcols = c("households", "individuals"))
#' y2 <- aggregate_duplicates(x, sumcols = c("households", "individuals"))
#' y3 <- aggregate_duplicates(x, idcols = c("municipality", "x", "y"))
#' y1
#' identical(y1, y2)
#' identical(y2, y3)
#' # Additional columns are dropped if not (possibly implicitly) referred to:
#' xx <- cbind(x, junk = 1:3)
#' xx
#' aggregate_duplicates(xx, idcols = c("municipality", "x", "y"),
#'                      sumcols = c("households", "individuals"))
aggregate_duplicates <- function(x, idcols = NULL, sumcols = NULL){
  stopifnot(requireNamespace("data.table"))
  if(is.null(idcols) && is.null(sumcols)){
    stop("Supply at least one of idcols and sumcols")
  }
  nam <- names(x)
  if(is.null(idcols)){
    if(is.null(sumcols) || !is.character(sumcols)){
      stop("When idcols is NULL sumcols must be a non-NULL character vector")
    }
    idcols <- nam[!(nam %in% sumcols)]
  }
  if(is.null(sumcols)){
    if(is.null(idcols) || !is.character(idcols)){
      stop("When sumcols is NULL idcols must be a non-NULL character vector")
    }
    sumcols <- nam[!(nam %in% idcols)]
  }
  nam <- nam[nam %in% c(idcols, sumcols)]
  x$id <- do.call(paste, append(as.list(x[, idcols]), list(sep = ",")))
  x <- data.table::as.data.table(x)
  x <- x[, lapply(.SD, sum), by = id, .SDcols = sumcols]
  x <- x[, paste(idcols) := data.table::tstrsplit(id, ",", fixed=TRUE)]
  x <- x[, id:=NULL]
  data.table::setcolorder(x, nam)
  for(col in c("x", "y")){
    data.table::set(x, j=col, value = as.numeric(x[[col]]))
  }
  return(as.data.frame(x))
}
