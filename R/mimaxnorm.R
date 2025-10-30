#' Min-Max Normalization Function
#' @param x A numeric vector to be normalized.
#' @return A numeric vector with values normalized to the range [0, 1].
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' minmaxnorm(x)
#' @export
minmaxnorm=function(x){(x-min(x))/(max(x)-min(x))}
check_nz=function(m){
  sapply(m,function(col){length(nonzero(as.numeric(col)))}) %>% as.data.frame() %>% dplyr::arrange(-.)
}