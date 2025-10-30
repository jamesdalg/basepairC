#' Convert to Normalized Counts
#'
#' This function adjusts the counts in a dataset by normalizing them based on the effective library size (ELS) and relative ELS.
#' It aims to reduce the impact of samples with high library size to give a fair comparison between condtions with different library sizes.
#'
#' @param dgel A data frame or list containing the fields `samples`, `norm.factors`, `lib.size`, and `counts`.
#'             `samples` should be a list or data frame that includes `norm.factors` and `lib.size` for each sample.
#'             `counts` should be a matrix or data frame of raw counts.
#' @return A matrix of normalized counts.
#' @examples
#' # Assuming `dgel` is your dataset with the appropriate structure:
#' normalized_counts <- convert_to_normalized_counts(dgel)
#' @references
#' For more information on the normalization process and its rationale, see the discussion at:
#' https://support.bioconductor.org/p/132240/
#' @export
convert_to_normalized_counts <- function(dgel) {
  # citation: https://support.bioconductor.org/p/132240/
  #another citation: https://support.bioconductor.org/p/124180/
  #yet another citation: https://support.bioconductor.org/p/73844/
  #and another: https://support.bioconductor.org/p/87185/
  els <- dgel$samples$norm.factors * dgel$samples$lib.size # get the ELS
  rels <- els / mean(els) # get the relative ELS
  norm <- t(t(dgel$counts) / rels) # divide by relative ELS to adjust the counts.
  # this will reduce the samples with high library size and increase those with low library size (below the mean).
  return(norm)
}
