#' Plot of the magnitude mapping
#'
#' @param x \code{'eqldensmap'} object
#' @param ... additional parameters to \code{\link{plot}}
#'
#' @export
#' @seealso \code{\link{eqldens}}
#'
#' @examples
#' M.ref <- runif(10,5,8)
#' plot(eqldens(1:10, M.ref))
plot.eqldensmap <- function(x, ...){
  data <- x[['mapped']]
  plot(data, ...)
}

#' Returns empirical cumulative distribution mappings
#'
#' @param x object
#' @param at numeric; the magnitudes to return the mapping for. If \code{NULL} (the default)
#' the function uses \code{\link{pretty}} to generate these.
#' @param inverse logical; should \code{at} be values of the ecdf? If \code{TRUE}
#' the ecdf is inverted to solve for the magnitudes
#' @param ... additional parameters to \code{\link{pretty}}
#'
#' @return a \code{data.frame} with the scale and mapping information
#' @export
#'
#' @examples
#' M.ref <- runif(10,5,8)
#' M.map <- eqldens(1:10, M.ref)
#' getScale(M.map)
getScale <- function(x, ...) UseMethod('getScale')

#' @rdname getScale
#' @export
getScale.eqldensmap <- function(x, at=NULL, inverse=FALSE, ...){
  Ref <- x[['reference']]
  CF <- Ref[['CF']]
  no.at <- is.null(at)
  if (no.at){
    at <- pretty(Ref[['data']], ...)
  } else {
    at <- na.omit(at)
    .inv_ecdf <- function(f){
      x <- environment(f)$x
      y <- environment(f)$y
      stats::approxfun(y, x)
    }
    if (inverse){
      CF <- .inv_ecdf(CF)
    }
  }
  values <- CF(at)
  map <- if (inverse & !no.at){
    data.frame(mapping=values, at=at)
  } else {
    data.frame(mapping=at, at=values)
  }
  return(map)
}
