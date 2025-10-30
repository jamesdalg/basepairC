#' Normalize a count matrix by removing zero rows/columns and applying SCN normalization
#'
#' This function processes a given count matrix by removing zero rows and columns,
#' applying SCN normalization, and then restoring the matrix to its original dimensions
#' and scale. Zero rows and columns are added back after normalization to maintain
#' the original structure of the count matrix.
#'
#' @param count_mat A matrix containing the counts to be normalized.
#' @param original_count_mat A matrix containing the original counts, used to revert the scaling after SCN normalization.
#' @return A normalized count matrix with the same dimensions and row/column names as the original matrix.
#' @importFrom Rfast rowsums colsums
#' @importFrom HiCcompare SCN
#' @examples
#' count_mat <- matrix(sample(0:10, 100, replace = TRUE), nrow = 10)
#' original_count_mat <- count_mat
#' result <- normalize_count_matrix(count_mat, original_count_mat)
#' @export
scn_normalize_to_original_scale=function(count_mat,original_scale=T){
original_count_mat <- count_mat
message_parallel("finding zero rows")
which(Rfast::rowsums(count_mat) == 0) -> zero_rows
which(Rfast::colsums(count_mat) == 0) -> zero_cols
basepairC::message_parallel("setting all_row_nums")
1:nrow(count_mat) -> all_row_nums
1:ncol(count_mat) -> all_col_nums
if (length(zero_rows) == 0) {
  non_zero_rows <- all_row_nums
} else {
  all_row_nums[-zero_rows] -> non_zero_rows
}
if (length(zero_cols) == 0) {
  non_zero_cols <- all_col_nums
} else {
  all_col_nums[-zero_cols] -> non_zero_cols
}
basepairC::message_parallel("subsetting count mat to remove zero rows and columns")
# all_row_nums[-zero_rows]->non_zero_rows
# all_col_nums[-zero_cols]->non_zero_cols
count_mat_without_zeros <- count_mat[non_zero_rows, non_zero_cols]
# count mat dim:
basepairC::message_parallel(paste0("count mat dim ", dim(count_mat)))
basepairC::message_parallel(paste0("count mat without zeros dim ", dim(count_mat_without_zeros)))
# run scn
scn_count_mat <- HiCcompare::SCN(count_mat_without_zeros + 1)
basepairC::message_parallel(paste0("scn count mat dim ", dim(scn_count_mat)))
# revert to original scale:
if(original_scale){
scn_count_mat <- minmaxnorm(scn_count_mat) * (max(original_count_mat) - min(nonzero(original_count_mat))) + min(nonzero(original_count_mat))}
#if the min is less than zero, add the min to all values to make them positive
if (min(nonzero(original_count_mat)) < 0) {
  scn_count_mat <- scn_count_mat + abs(min(nonzero(original_count_mat)))
}
# add the zero rows back in, in the correct order.
#    rbind(scn_count_mat,matrix(0,nrow = length(zero_rows),ncol=ncol(scn_count_mat)) ) -> scn_count_mat_with_zeros
# add the zero columns back in, in the correct order.
#    cbind(scn_count_mat_with_zeros,matrix(0,nrow = nrow(scn_count_mat_with_zeros),ncol=length(zero_cols)) ) -> scn_count_mat_with_zeros
# reorder the rows and columns to their original order.
# create an empty matrix of the same dimensions as the original count mat
matrix(0, nrow = nrow(count_mat), ncol = ncol(count_mat)) -> scn_count_mat_with_zeros
basepairC::message_parallel(paste0("scn_count_mat_with_zeros dim ", dim(scn_count_mat_with_zeros)))
# fill in the non-zero rows and columns with the scn count mat
scn_count_mat_with_zeros[non_zero_rows, non_zero_cols] <- scn_count_mat
count_mat <- scn_count_mat_with_zeros
basepairC::message_parallel(paste0("scn_count_mat_with_zeros dim ", dim(scn_count_mat_with_zeros)))
basepairC::message_parallel(paste0("final count_mat dim ", dim(count_mat)))
basepairC::message_parallel(paste0("length of rownames of original_count_mat ", rownames(original_count_mat)))
basepairC::message_parallel(paste0("length of colnames of original_count_mat ", colnames(original_count_mat)))
rownames(count_mat) <- rownames(original_count_mat)
colnames(count_mat) <- colnames(original_count_mat)
return(count_mat)
}