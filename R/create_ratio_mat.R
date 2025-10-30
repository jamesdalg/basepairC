#' Create Ratio Matrix Function
#' 
#' This function takes two matrices, performs various operations on them, and
#' returns a ratio matrix. It includes options for capping values, applying a
#' threshold, zeroing the top values, and adjusting for total reads.
#'
#' @param mat1 A numeric matrix.
#' @param mat2 A numeric matrix.
#' @param cap The cap value for thresholding (default is 100).
#' @param thresh The threshold value for thresholding.
#' @param zero_top Logical, indicating whether to zero out top values (default is FALSE).
#' @param quantile The quantile value for thresholding (default is 0).
#' @param adjust_for_total_reads A character string, specifying whether to adjust for total reads "before_cleaning" or "after_cleaning".
#'
#' @return A numeric matrix with ratio values.
#'
#' @examples
#' mat1 <- matrix(1:10, nrow = 2)
#' mat2 <- matrix(11:20, nrow = 2)
#' create_ratio_mat(mat1, mat2)
#'
#' @export
create_ratio_mat=function(mat1,mat2,cap=100,thresh,zero_top=F,quantile=0,adjust_for_total_reads="after_cleaning"){
  zero_top->zt
  ct_mat1=cap_thresh(mat1,cap = cap,thresh=thresh,zero_top=zt,quantile = quantile)
  ct_mat2=cap_thresh(mat2,cap = cap,thresh=thresh,zero_top=zt,quantile = quantile)
  total_read_ratio=sum(ct_mat1)/sum(ct_mat2)
  if(sum(ct_mat1)>=sum(ct_mat2)){ct_mat1=ct_mat1/total_read_ratio} else {ct_mat1=ct_mat1*total_read_ratio }
  ct_mat1/ct_mat2->ratio_mat
  ratio_mat[is.na(ratio_mat)]=0
  ratio_mat[is.infinite(ratio_mat)]=0
  log2(ratio_mat+1)->log2_ratio_mat
  return(log2_ratio_mat)
}
