#' Perform Read Correction on a List of Matrices
#'
#' Adjusts the values in each matrix of the list (except for the one with the smallest
#' total read count) to match the total sum of the matrix with the smallest read count.
#' The adjustment is based on the ratio of total sums between each matrix and the one
#' with the smallest total sum. This function aims to equalize the total sums of all
#' matrices with that of the matrix with the smallest total sum.
#'
#' @param matrices A list of numeric matrices of the same dimensions. Each matrix
#'        will be adjusted relative to the matrix with the smallest total sum.
#'
#' @return A list containing the adjusted matrices.
#'
#' @examples
#' mat1 <- matrix(c(1, 2, 3, 4), nrow=2)
#' mat2 <- matrix(c(2, 4, 6, 8), nrow=2)
#' mat3 <- matrix(c(3, 6, 9, 12), nrow=2)
#' adjusted_matrices <- total_read_correction_list(list(mat1, mat2, mat3))
#' adjusted_matrices[[1]] # Adjusted mat1
#' adjusted_matrices[[2]] # Adjusted mat2
#' adjusted_matrices[[3]] # Original mat3 (smallest read count matrix in this case)
total_read_correction_list <- function(matrices) {
  # Identify the matrix with the smallest total read count
  total_sums <- sapply(matrices, sum)
  smallest_matrix_idx <- which.min(total_sums)
  smallest_matrix <- matrices[[smallest_matrix_idx]]
  
  # Adjust each matrix relative to the smallest read count matrix
  adjusted_matrices <- lapply(matrices, function(mat) {
    if (identical(mat, smallest_matrix)) {
      return(mat) # Return the smallest matrix unmodified
    }
    total_read_ratio <- sum(mat) / sum(smallest_matrix)
    if (sum(mat) >= sum(smallest_matrix)) {
      return(mat / total_read_ratio)
    } else {
      return(mat * total_read_ratio)
    }
  })
  
  return(adjusted_matrices)
}