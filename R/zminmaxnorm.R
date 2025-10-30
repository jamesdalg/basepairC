#' Z-Score Min-Max Normalization Function
#' 
#' This function takes a numeric vector, standardizes it (z-score normalization),
#' and then performs min-max normalization on the standardized values, scaling them
#' to the range [0, 1].
#'
#' @param x A numeric vector to be normalized.
#'
#' @return A numeric vector with z-score min-max normalized values.
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' zminmaxnorm(x)
#'
#' @export
zminmaxnorm=function(x){minmaxnorm((x-mean(x))/sd(x))}


#' Finite Elements Function
#'
#' This function takes a numeric vector and returns only the finite elements.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector with finite values.
#'
#' @examples
#' x <- c(1, 2, Inf, 3, NaN)
#' noninf(x)
#'
#' @export
noninf=function(x){x[is.finite(x)]}

#' Z-Transformation Function
#'
#' This function takes a numeric vector and performs z-transformation.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector with z-transformed values.
#'
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' ztrans(x)
#'
#' @export
ztrans=function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)}