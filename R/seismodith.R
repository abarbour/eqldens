
#' Dither earthquake magnitudes assuming a Gutenberg-Richter distribution
#'
#' @param M numeric; earthquake magnitude
#' @param dither.by numeric; the amoung to dither by
#' @param b.value numeric; the 'b-value' in the Gutenberg-Richter distribution
#' @param offset numeric; a constant to add to all magnitudes
#'
#' @return a \code{data.frame} with the original magnitudes, the dithered
#' magnitudes, and the uniform distribution used in the dithering
#' @export
#' @seealso \code{\link{eqldens}}
#'
#' @examples
#' 1
seismodith <- function(M, dither.by=0.1, b.value=1.0, offset=0.0){
	M <- Mo <- as.numeric(M)
	nm <- length(M)
	offset <- offset - dither.by/2	# offs=offs-del/2
	b.value <- b.value * log(10)		# b=b*alog(10.)
	F. <- 1 - exp(-b.value * dither.by)	# f=1-exp(-b*del)
	u <- runif(nm, min=0, max=1)		# [ read-loop: rand() ]
	M <- M - log(1 - u*F.)/b.value + offset	# rmag <- rmag-alog(1-u*f)/b+offs
	data.frame(Mo=Mo, Mq=M, U=u)
}

#' Map magnitudes to a reference distribution
#'
#' @param M numeric; the magnitudes to map
#' @param Mref numeric; the reference magnitudes used to create the mapping
#'
#' @return An object with class \code{'eqldensmap'}
#' @export
#' @seealso \code{\link{seismodith}}
#'
#' @examples
#' M.ref <- runif(10,5,8)
#' eqldens(1:10, M.ref)
eqldens <- function(M, Mref){

	M <- as.numeric(M)
	Mref <- as.numeric(M)

	# Create a function to return empirical cdf
	CF <- stats::ecdf(Mref)
	#kts <- stats::knots(CF)
	# and then a monotone cubic spline at the cdf
	KCF <- stats::splinefun(CF(Mref), method="monoH.FC")

	# Now 'map' the magnitudes onto this reference distro
	Mmap <- CF(M)

	# and return
	Results <- list(
		reference = list(
			data = Mref,
			CF = CF,
			KCF = KCF
		),
		mapped = data.frame(M, Mmap)
	)
	class(Results) <- c('eqldensmap', 'list')
	return(Results)
}
