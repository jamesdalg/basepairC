#' Correct Mappability and GC Bias in Count Matrices
#'
#' Use a negative binomial generalized linear model to adjust counts based on mappability and GC content.
#'
#' @param original_mat A matrix of original counts to be corrected. Each row and column should represent different genomic features or samples.
#' @param dt_with_map_and_gc A data table containing mappability and GC content information for genomic regions. The table should have columns for chromosome (`chr`), start and end positions (`start`, `end`), mappability (`map_24bp`), and GC content (`gc_1bp`).
#'
#' @return A matrix of counts corrected for mappability and GC content biases, with the same dimensions and row/column names as the input matrix.
#'
#' @examples
#' # Assuming `original_mat` is your matrix of counts and `dt_with_map_and_gc` is your data table with mappability and GC content:
#' corrected_mat <- correct_mappability_gc_bias(original_mat, dt_with_map_and_gc)
#'
#' @references
#' Hu M, Deng K, Selvaraj S, Qin ZS, Ren B and Liu JS (2012) HiCNorm: removing biases in Hi-C data via Poisson regression. (2012) Bioinformatics 28 (23), 3131-3133.
#'
#' @export
correct_mappability_gc_bias <- function(original_mat, dt_with_map_and_gc) {
  # browser()
  print("melting matrix, splitting into cols, getting averages, and merging with gene_locus_dt_bp.")
  original_mat %>%
    reshape2::melt() %>%
    dtplyr::lazy_dt() %>%
    dplyr::rename(loc1 = Var1, loc2 = Var2, counts = value) %>%
    tidyr::separate(loc1, sep = "_", into = c("chr1", "start1", "end1")) %>%
    tidyr::separate(loc2, sep = "_", into = c("chr2", "start2", "end2")) %>%
    dplyr::mutate(start1 = as.numeric(start1), start2 = as.numeric(start2), end1 = as.numeric(end1), end2 = as.numeric(end2)) %>%
    data.table::as.data.table() %>%
    data.table::merge.data.table(dt_with_map_and_gc, by.x = c("chr1", "start1", "end1"), by.y = c("chr", "start", "end")) %>%
    data.table::merge.data.table(dt_with_map_and_gc, by.x = c("chr2", "start2", "end2"), by.y = c("chr", "start", "end")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(mappability = log(map_24bp.x * map_24bp.y), gc = log(gc_1bp.x * gc_1bp.y)) -> mat_norm_filter_melted
  MASS::glm.nb(data = mat_norm_filter_melted, formula = counts ~ gc + offset(mappability)) -> mat_norm_filter_model
  round(mat_norm_filter_melted$counts / exp(mat_norm_filter_model$coefficients[1] + mat_norm_filter_model$coefficients[2] * mat_norm_filter_melted$gc + mat_norm_filter_melted$mappability)) -> mat_norm_filter_melted$gc_map_corr_counts
  mat_norm_filter_melted$gc_map_corr_counts %>%
    as.numeric() %>%
    matrix(nrow = nrow(original_mat), ncol = ncol(original_mat)) -> ctrl_mat_norm_filter_cleaned_counts_mat
  colnames(ctrl_mat_norm_filter_cleaned_counts_mat) <- colnames(original_mat)
  rownames(ctrl_mat_norm_filter_cleaned_counts_mat) <- rownames(original_mat)
  return(ctrl_mat_norm_filter_cleaned_counts_mat)
}
