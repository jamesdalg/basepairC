#' Perform SCN and Min-Max Normalization on a Count Matrix
#'
#' This function applies symmetric normalization (SCN) and min-max normalization to a given count matrix.
#' It first makes the matrix symmetric if it is square, removes rows and columns with all zeros, applies SCN,
#' then applies min-max normalization to scale the values back to the original scale. It finally reintegrates
#' the zero rows and columns in their original positions and returns the normalized matrix with original row
#' and column names.
#'
#' @param count_mat A numeric matrix representing the count data to be normalized.
#'
#' @return A numeric matrix with the same dimensions as `count_mat`, where the data has been symmetrically
#' normalized and scaled using min-max normalization, with zero rows and columns preserved in their original
#' order.
#'
#' @examples
#' # Generate a random count matrix
#' set.seed(123)
#' count_mat <- matrix(rnorm(100), nrow = 10)
#' # Apply SCN and min-max normalization
#' normalized_mat <- scn_minmaxnorm(count_mat)
#' @references Performs Sequential Component Normalization as described by Cournac. Coded using details in the manuscript. Cournac A, Marie-Nelly H, Marbouty M, Koszul R, Mozziconacci J. Normalization of a chromosomal contact map. BMC Genomics. 2012;13: 436. doi:10.1186/1471-2164-13-436
#' John C. Stansfield, Kellen G. Cresswell, Vladimir I. Vladimirov, Mikhail G. Dozmorov, HiCcompare: an R-package for joint normalization and comparison of HI-C datasets. BMC Bioinformatics. 2018 Jul 31;19(1):279. doi: 10.1186/s12859-018-2288-x.
#'
#'
#' @importFrom HiCcompare SCN
#' @importFrom Rfast rowsums colsums
#' @export
scn_minmaxnorm <- function(count_mat) {
  # browser()
  original_count_mat <- count_mat
  if (nrow(count_mat) == ncol(count_mat)) {
    count_mat %>% basepairC:::make_symmetric() -> count_mat
  }

  which(Rfast::rowsums(count_mat) == 0) -> zero_rows
  which(Rfast::colsums(count_mat) == 0) -> zero_cols
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
  count_mat_without_zeros <- count_mat[non_zero_rows, non_zero_cols]
  # run scn
  scn_count_mat <- HiCcompare::SCN(count_mat_without_zeros)
  # revert to original scale:

  scn_count_mat <- minmaxnorm(scn_count_mat) * (max(original_count_mat) - min(nonzero(original_count_mat))) + min(nonzero(original_count_mat))
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
  # fill in the non-zero rows and columns with the scn count mat
  scn_count_mat_with_zeros[non_zero_rows, non_zero_cols] <- scn_count_mat
  count_mat <- scn_count_mat_with_zeros
  rownames(count_mat) <- rownames(original_count_mat)
  colnames(count_mat) <- colnames(original_count_mat)
  return(count_mat)
}
