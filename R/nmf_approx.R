#' Approximate a Matrix using Non-negative Matrix Factorization (NMF)
#'
#' This function approximates a given matrix using Non-negative Matrix Factorization (NMF).
#'
#' @param mat A numeric matrix. The input matrix to be approximated.
#' @param N An integer. The number of components to use for the NMF approximation.
#' @param nmf_method A character string. The method to use for the NMF algorithm. Default is "lee".
#'
#' @return A numeric matrix. The approximated matrix.
#'
#' @details This function performs Non-negative Matrix Factorization (NMF) on the input matrix `mat`
#' using the specified number of components `N` and the specified method `nmf_method`. It returns the
#' approximated matrix obtained by multiplying the basis matrix `W` and the coefficient matrix `H`.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(runif(100), nrow = 10)
#' N <- 3
#' approx_mat <- nmf_approx(mat, N)
#' }
#'
#' @importFrom NMF nmf basis coef
#' @export
nmf_approx <- function(mat, N,nmf_method="lee"){
  res <- NMF::nmf(mat, N, method = nmf_method)
  W <- NMF::basis(res)
  H <- NMF::coef(res)
  mat.approx <- W %*% H
  return(mat.approx)
}