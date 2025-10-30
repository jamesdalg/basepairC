
#' Non-Zero Elements Function
#'
#' This function takes a numeric vector and returns only the non-zero elements.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector with non-zero values.
#'
#' @examples
#' x <- c(0, 1, 2, 0, 3)
#' nonzero(x)
#'
#' @export
nonzero=function(x){x[x!=0]}