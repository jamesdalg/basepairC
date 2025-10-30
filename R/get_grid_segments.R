#' Get grid segments from a sample matrix list
#'
#' This function calculates grid breakpoints for the rows and columns of a matrix, based on a specified interval (`every_bp`). 
#' The breakpoints represent the positions at which segments of the matrix are defined. The function returns a list with 
#' breakpoints for rows, columns, and their transposed equivalents.
#'
#' @param sample_mat_list A list of matrices. The function assumes all matrices have the same dimensions, and it uses the first matrix in the list to determine the dimensions.
#' @param every_bp An integer specifying the interval (in base pairs) at which to define breakpoints for the matrix grid. Default is 500.
#'
#' @return A list containing four elements:
#' \describe{
#'   \item{breakpoints_col}{A vector of breakpoints for the columns of the matrix.}
#'   \item{breakpoints_row}{A vector of breakpoints for the rows of the matrix.}
#'   \item{t_breakpoints_col}{A vector of breakpoints for the columns of the transposed matrix.}
#'   \item{t_breakpoints_row}{A vector of breakpoints for the rows of the transposed matrix}
#' }
#'
#' @examples
#' # Example usage:
#' mat1 <- matrix(1:10000, nrow = 100, ncol = 100)
#' mat_list <- list(mat1)
#' seg_obj <- get_grid_segments(mat_list, every_bp = 20)
#' print(seg_obj)
#'
#' @export
get_grid_segments=function(sample_mat_list,every_bp=500){
  dim(sample_mat_list[[1]])->dim_mat
  seq(1,dim_mat[1],every_bp)->grid_seg_num_col
  setdiff(c(dim_mat[1],grid_seg_num_col),dim_mat[1])->grid_seg_num_col
  seq(1,dim_mat[2],every_bp)->grid_seg_num_row
  setdiff(c(dim_mat[2],grid_seg_num_row),dim_mat[2])->grid_seg_num_row
  seg_obj=list(breakpoints_col=grid_seg_num_col,breakpoints_row=grid_seg_num_row,t_breakpoints_col=grid_seg_num_row,t_breakpoints_row=grid_seg_num_col)
  return(seg_obj)
}