#' Fused Lasso Filtering Function
#' 
#' This function applies the fused lasso filtering technique to a given matrix.
#'
#' @param mat A numeric matrix to be filtered.
#' @param lambda A numeric value representing the regularization parameter.
#' @param eps Convergence threshold paramter (see glmgen::fusedLattice).
#' @param maxIter  maximum number of iterations for convergence.
#'
#' @return A filtered numeric matrix.
#'
#' @examples
#' mat <- matrix(1:10, nrow = 2)
#' fused_lasso_filter(mat, lambda = 0.1, eps = 0.001, maxIter = 100)
#'
#' @export
fused_lasso_filter=function(mat,lambda,eps,maxIter){
  glmgen::fusedLattice(mat,lambda=lambda,eps=eps,maxIter=maxIter)->z
  t(mat)%*%matrix(unlist(z$beta),ncol=ncol(mat))->yhat_matrix
  return(yhat_matrix)
}
