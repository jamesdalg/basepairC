#' Perform 2D LASSO regularization on a matrix
#'
#' This function applies 2D LASSO regularization on a given matrix using the fusedLattice function from the glmgen package. It returns the regularized matrix and the beta matrix used for regularization. Optionally, it can remove elements of the matrix that are zero after regularization.
#'
#' @param mat A numeric matrix on which 2D LASSO regularization is to be applied.
#' @param lambda Regularization parameter for LASSO. Defaults to 5.
#' @param eps Convergence threshold for the algorithm. Defaults to 1e-2.
#' @param maxIter Maximum number of iterations for the algorithm. Defaults to 20.
#' @param remove_zero Logical indicating whether to remove zero elements from the regularized matrix. Defaults to FALSE.
#'
#' @return A list containing the regularized matrix (`lasso_mat`) and the beta matrix (`beta_mat`). If `remove_zero` is TRUE, `lasso_mat` will only contain non-zero elements.
#'
#' @examples
#' mat <- matrix(rnorm(100), 10, 10)
#' result <- lasso2d(mat)
#' str(result)
#'
#' @export
lasso2d=function(mat,lambda=5,eps=1e-2,maxIter=20,remove_zero=F){
  beta_mat=matrix(glmgen::fusedLattice(as.matrix(mat), lambda = lambda, eps=eps, maxIter=maxIter)$beta,nrow=nrow(mat))
  t_beta_mat=matrix(glmgen::fusedLattice(as.matrix(t(mat)), lambda = lambda, eps=eps, maxIter=maxIter)$beta,nrow=ncol(mat))
  as.matrix(beta_mat*mat*t(t_beta_mat))->lasso_mat
  if(remove_zero){
    return(list(lasso_mat=lasso_mat[lasso_mat>0],beta_mat=beta_mat))
  }
  return(list(lasso_mat=lasso_mat,beta_mat=beta_mat))
}
sumlogp=function(p){
  as.numeric(metap::sumlog(p)$p)
}