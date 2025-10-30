#' Get Sum of Segments from Sample Matrix List
#'
#' This function computes the sum of matrices from a list of sample matrices and processes
#' the resulting matrix by applying a segmentation preprocessing step.
#'
#' @param sample_mat_list A list of matrices, where each matrix corresponds to a sample.
#' @param n_segs An integer representing the number of segments for both rows and columns in the Lasso segmentation (default is 9).
#'
#' @return A list containing:
#' \describe{
#'   \item{lasso_segments}{Lasso-based segmentation results.}
#'   \item{orig_segments}{Original segmentation results.}
#' }
#' @examples
#' \dontrun{
#'   get_sum_segments(sample_mat_list, n_segs = 9)
#' }
#' @export
get_sum_segments=function(sample_mat_list,n_segs=9){
  Reduce("+",sample_mat_list)->sum_matrix_all_samples
  preprocess_mat(sum_matrix_all_samples,verbose=T,dsfactor=1,lasso_seg_max_col = n_segs,lasso_seg_max_row = n_segs)->sum_matrix_all_samples_df
  return(list(lasso_segments=attr(sum_matrix_all_samples_df,"segments_lasso_obj"),orig_segments=attr(sum_matrix_all_samples_df,"segments_orig_obj")))
}