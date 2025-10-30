#' Perform Read Correction on a Pair of Matrices
#'
#' This function performs read correction on a pair of matrices by adjusting the
#' values in the first matrix (`mat1`) to match the total sum of the second matrix
#' (`mat2`). The adjustment is based on the ratio of total reads in `mat1` and
#' `mat2`. The function aims to equalize the total sum of `mat1` with that of
#' `mat2` by either dividing or multiplying its elements by the calculated total
#' read ratio.
#'
#' @param mat1 A numeric matrix. The first matrix to be adjusted.
#' @param mat2 A numeric matrix. The second matrix whose total sum is used as the
#'             target for adjustment of `mat1`.
#'
#' @return A list containing two elements: the adjusted `mat1` and the original
#'         `mat2`, in that order.
#'         
#' @examples
#' mat1 <- matrix(c(1, 2, 3, 4), nrow=2)
#' mat2 <- matrix(c(2, 4, 6, 8), nrow=2)
#' corrected_matrices <- total_read_correction(mat1, mat2)
#' corrected_matrices[[1]] # Adjusted mat1
#' corrected_matrices[[2]] # Original mat2
total_read_correction=function(mat1,mat2) {
  total_read_ratio=sum(mat1)/sum(mat2)
  if(sum(mat1)>=sum(mat2)){mat1=mat1/total_read_ratio} else {mat1=mat1*total_read_ratio }
  return(list(mat1,mat2))
}
