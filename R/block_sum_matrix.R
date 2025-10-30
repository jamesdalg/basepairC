#' Downsample a matrix by summing blocks of values
#' 
#' @param matrix The input matrix to downsample
#' @param downsample_factor_rows The downsampling factor for rows
#' @param downsample_factor_cols The downsampling factor for columns (default: \code{downsample_factor_rows})
#' @param average perform average instead of sum
#' @return A downsampled matrix
#' 
#' @details This function downsample a matrix by summing blocks of values. It takes an input matrix and two downsampling factors for rows and columns. It then creates downsampling matrices for rows and columns based on the provided factors. It then performs the downsampling operation by multiplying the input matrix with the transpose of the row downsampling matrix, then with the input matrix again, and finally with the column downsampling matrix. The result is a downsampled matrix.
#' 
#' @examples
#' matrix(1, nrow = 5, ncol = 4) %>% block_sum_matrix(5, 2)
#' 
#' @export
block_sum_matrix <- function(matrix, downsample_factor_rows, downsample_factor_cols = downsample_factor_rows,average=F) {
  # Check that the dimensions of the matrix are divisible by the downsampling factors
  if (nrow(matrix) %% downsample_factor_rows != 0 || ncol(matrix) %% downsample_factor_cols != 0) {
    stop("The dimensions of the matrix must be divisible by the downsampling factors")
  }
  
  # Create grouping variables for the rows and columns
  row_groups <- rep(seq_len(nrow(matrix) / downsample_factor_rows), each = downsample_factor_rows)
  col_groups <- rep(seq_len(ncol(matrix) / downsample_factor_cols), each = downsample_factor_cols)
  
  # Sum the rows of the matrix based on the row groups
  #browser()
  row_sums <- rowsum(as.matrix(matrix), row_groups)
  
  # Transpose the matrix, sum the columns (which are the original rows) based on the column groups, and transpose the result
  downsampled_matrix <- t(rowsum(t(row_sums), col_groups))
  if(average){downsampled_matrix=downsampled_matrix/(downsample_factor_rows*downsample_factor_cols)}
  return(downsampled_matrix)
}
# block_sum_matrix <- function(matrix, downsample_factor_rows, downsample_factor_cols = downsample_factor_rows) {
#   # Check that the dimensions of the matrix are divisible by the downsampling factors
#   if (nrow(matrix) %% downsample_factor_rows != 0 || ncol(matrix) %% downsample_factor_cols != 0) {
#     stop("The dimensions of the matrix must be divisible by the downsampling factors")
#   }
#   
#   # Calculate the dimensions of the downsampled matrix
#   num_rows <- nrow(matrix) / downsample_factor_rows
#   num_cols <- ncol(matrix) / downsample_factor_cols
#   
#   # Reshape the matrix into a 3D array and sum over the third dimension
#   array_3d <- array(matrix, dim = c(downsample_factor_rows, downsample_factor_cols, num_rows * num_cols))
#   downsampled_matrix <- Rfast::colsums(apply(array_3d, 3, sum))
#   
#   # Reshape the result into a matrix
#   downsampled_matrix <- matrix(downsampled_matrix, nrow = num_rows, ncol = num_cols)
#   
#   return(downsampled_matrix)
# }
# block_sum_matrix <- function(matrix, downsample_factor_rows, downsample_factor_cols = downsample_factor_rows) {
#   # Check that the dimensions of the matrix are divisible by the downsampling factors
#   if (nrow(matrix) %% downsample_factor_rows != 0 || ncol(matrix) %% downsample_factor_cols != 0) {
#     stop("The dimensions of the matrix must be divisible by the downsampling factors")
#   }
# 
#   # Create grouping variables for the rows and columns
#   row_groups <- rep(seq_len(nrow(matrix) / downsample_factor_rows), each = downsample_factor_rows)
#   col_groups <- rep(seq_len(ncol(matrix) / downsample_factor_cols), each = downsample_factor_cols)
# 
#   # Sum the rows of the matrix based on the row groups
#   row_sums <- Rfast::rowsums(matrix, row_groups)
# 
#   # Transpose the matrix, sum the columns (which are the original rows) based on the column groups, and transpose the result
#   downsampled_matrix <- t(Rfast::rowsums(t(row_sums), col_groups))
# 
#   return(downsampled_matrix)
# }
# block_sum_matrix <- function(matrix, downsample_factor_rows, downsample_factor_cols = downsample_factor_rows) {
#   # Create the downsampling matrices for rows and columns
#   create_downsample_matrix <- function(n, r) {
#     kronecker(diag(n / r), matrix(1, nrow = r, ncol = r))
#   }
#   
#   # Row downsampling matrix and its transpose
#   row_mat <- create_downsample_matrix(nrow(matrix), downsample_factor_rows)
#   t_row_mat <- t(row_mat)
#   
#   # Column downsampling matrix (use row matrix if the factor is the same)
#   if (missing(downsample_factor_cols) || downsample_factor_cols == downsample_factor_rows) {
#     col_mat <- row_mat
#   } else {
#     col_mat <- create_downsample_matrix(ncol(matrix), downsample_factor_cols)
#   }
#   t_col_mat <- t(col_mat)
#   
#   # Perform the downsampling
#   browser()
#   downsampled_matrix <- t_row_mat %*% matrix %*% t_col_mat
#   
#   return(downsampled_matrix)
# }
downsample_with_names=function(mat,factor,mean=F){
  basepairC::block_sum_matrix(mat,downsample_factor_rows = factor,downsample_factor_cols = factor,average=mean)->block_10bp
  CNVScope::downsample_genomic_matrix(mat,downsamplefactor = factor)->ds_10bp
  colnames(block_10bp)=colnames(ds_10bp)
  rownames(block_10bp)=rownames(ds_10bp)
  return(block_10bp)
}