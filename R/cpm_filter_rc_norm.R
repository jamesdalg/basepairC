#' CPM Filter Read Counts Normalization
#'
#' This function normalizes and filters read counts using a DESeq2-like method, correcting for conditions 
#' and other factors provided in the input object. It returns the filtered and normalized read counts.
#' *Input is the output from read_replicates()*
#'
#' @param rr A list containing the following components:
#' \describe{
#'   \item{sample_mat}{Matrix of read counts, where rows represent features (e.g., genomic regions) and columns represent samples.}
#'   \item{conditions}{A vector of conditions corresponding to the samples.}
#'   \item{loc1}{Numeric vector specifying the first set of genomic locations.}
#'   \item{loc2}{Numeric vector specifying the second set of genomic locations.}
#'   \item{rep_dim}{A vector of length 2 specifying the output matrix dimensions (number of rows and columns).}
#'   \item{underscored_positions_col}{A logical indicating whether to underscore column positions.}
#'   \item{underscored_positions_row}{A logical indicating whether to underscore row positions.}
#' }
#'
#' @return A list containing the normalized and filtered read count matrix (`rr_obj`).
#'
#' @details
#' The function calls `normalizeFilterReadCounts` to normalize read counts using the DESeq2 method 
#' without estimating dispersion, applying filtering based on differences, and optionally skipping matrices.
#' 
#' @seealso \code{\link{normalizeFilterReadCounts}}
#'
#' @examples
#' \dontrun{
#' rr <- list(
#'   sample_mat = matrix(rnorm(100), nrow=10, ncol=10),
#'   conditions = factor(rep(c("A", "B"), each=5)),
#'   loc1 = 1:10,
#'   loc2 = 11:20,
#'   rep_dim = c(10, 10),
#'   underscored_positions_col = TRUE,
#'   underscored_positions_row = FALSE
#' )
#' cpm_filter_rc_norm(rr)
#' }
#' 
#' @export
cpm_filter_rc_norm=function(rr){
  normalizeFilterReadCounts(sample_mat = rr$sample_mat,
                            conditions = rr$conditions,
                            loc1 = rr$loc1,
                            loc2 = rr$loc2,
                            nrow_output = rr$rep_dim[1],
                            ncol_output = rr$rep_dim[2],
                            output_type = "corrected_counts",
                            filter = T,
                            underscored_positions_col = rr$underscored_positions_col,
                            underscored_positions_row = rr$underscored_positions_row,
                            norm_factor_type = "DESeq2",
                            filter_type = "diff",
                            skip_mats = T,
                            estimate_dispersion = F,
                            norm_offsets = F)->normalizeFilterReadCounts_output_counts_no_scn_deseq2_diff_ds
  return(normalizeFilterReadCounts_output_counts_no_scn_deseq2_diff_ds$rr_obj)
}