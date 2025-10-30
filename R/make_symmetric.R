#' Enforce matrix symmetry
#' 
#' Adds the upper and lower triangular parts, after forcing symmetry on each
#' @keywords matrix symmetry symmetric upper lower
#' @importFrom Matrix forceSymmetric
#' @param mat A matrix
#' @return A symmetric matrix
make_symmetric=function(mat) {
  return(as.matrix(forceSymmetric(as.matrix(mat),"U"))+as.matrix(forceSymmetric(as.matrix(mat),"L")))
}