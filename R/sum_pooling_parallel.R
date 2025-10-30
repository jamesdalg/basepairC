#' Sum Pooling with Parallel Processing Function
#'
#' This function performs sum pooling on a matrix with parallel processing support.
#'
#' @param mat A numeric matrix to be pooled.
#' @param factor The pooling factor for aggregating matrix elements.
#' @param ncores The number of CPU cores to use for parallel processing (default is the number of detected cores).
#'
#' @return A pooled numeric matrix using sum pooling.
#'
#' @examples
#' mat <- matrix(1:16, nrow = 4)
#' sum_pooling_parallel(mat, factor = 2)
#'
#' @export
sum_pooling_parallel <- function(mat, factor,ncores=parallel::detectCores()) {
  if(factor==1){return(as.matrix(mat))}
  indices <- expand.grid(seq(1, nrow(mat), by=factor), seq(1, ncol(mat), by=factor))
  pooling_function <- function(idx) {
    i <- as.numeric(idx[1])
    j <- as.numeric(idx[2])
    return(sum(mat[i:(i+factor-1), j:(j+factor-1)]))
  }
  if(ncores==1){pooled_values <- pbapply::pblapply(1:nrow(indices), function(k) pooling_function(indices[k,]))} else{
    pooled_values <- pbmcapply::pbmclapply(1:nrow(indices), function(k) pooling_function(indices[k,]), mc.cores=ncores)}
  pooled_matrix <- matrix(unlist(pooled_values), nrow=sqrt(length(pooled_values)))
  return(pooled_matrix)
}