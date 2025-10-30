#' Equalize Sample Coverage
#'
#' This function takes a matrix of samples and equalizes the coverage by sampling non-zero entries.
#' It handles numeric and non-numeric columns separately, ensuring that the numeric columns have equalized
#' coverage across samples. The non-numeric columns are returned unchanged.
#'
#' @param sample_mat A data frame or matrix where columns represent samples. Columns can be a mix of numeric and non-numeric data.
#' @param n_samples An integer indicating the number of samples to generate in the equalization process.
#'
#' @return A data frame with the same dimensions as `sample_mat`, where numeric columns have had their coverage equalized
#'   across samples, and non-numeric columns are returned unchanged.
#'
#' @details 
#' The function first separates numeric and non-numeric columns. It then processes the numeric columns to equalize coverage 
#' by sampling non-zero entries in each column such that the number of non-zero entries matches the minimum across all columns.
#' Finally, it recombines the numeric and non-numeric columns and returns the result.
#'
#' @importFrom Rfast colsums rowmeans
#' @importFrom parallel mclapply detectCores
#' @importFrom data.table data.table
#' @importFrom dplyr bind_cols
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' sample_matrix <- data.frame(
#'   sample1 = c(1, 2, 0, 4),
#'   sample2 = c(0, 5, 6, 0),
#'   sample3 = c(7, 0, 9, 10),
#'   non_numeric = c("A", "B", "C", "D")
#' )
#' n_samples <- 100
#' equalized_samples <- equalize_sample_coverage(sample_matrix, n_samples)
#' }
#'
#' @export
equalize_sample_coverage=function(sample_mat,n_samples,max_cores=parallel::detectCores()){
  message_parallel("equalize_sample_coverage")
  sapply(sample_mat,is.numeric)->numeric_columns
  sample_mat=data.table::as.data.table(sample_mat)
  sample_mat[,!..numeric_columns]->sample_mat_nonnumeric
  sample_mat[,..numeric_columns]->sample_mat_numeric
  total_counts=Rfast::colsums(as.matrix(sample_mat_numeric),parallel=0)
  sample_mat_nz=sample_mat_numeric
  sample_mat_nz[sample_mat_nz!=0]<-1
  total_nz_counts=Rfast::colsums(as.matrix(sample_mat_nz),parallel=0)
  min(total_nz_counts)->min_total_nz_counts
  
  #options(mc.cores=min(n_samples,parallel::detectCores()))
  #max_cores=parallel::detectCores()
  outer_cores=min(ncol(sample_mat_numeric),max_cores)
  inner_cores=max(max_cores-outer_cores*n_samples,1)
    #min(n_samples,parallel::detectCores())
  #browser()
  parallel::mclapply(1:ncol(sample_mat_numeric),function(i){
    parallel::mclapply(1:n_samples,function(j){
      message_parallel(paste0("resample: ",j," Column: ",i))
      sample_col=as.matrix(sample_mat_numeric)[,i]
      nz_positions=which(sample_col!=0)
      sampled_positions=sample(nz_positions,min_total_nz_counts)
      unsampled_positions=setdiff(nz_positions,sampled_positions)
      sample_col[unsampled_positions]<-0
      return(sample_col)
    },mc.cores=inner_cores) %>% do.call(cbind,.) %>% Rfast::rowmeans()  -> equalized_sample_mat
    data.table::data.table(equalized_sample_mat)->equalized_sample_dt
    colnames(equalized_sample_dt)=colnames(sample_mat_numeric)[i]
    return(equalized_sample_dt)
  },mc.cores=outer_cores) %>% do.call(cbind,.) -> equalized_sample_dt_complete
  #add back the non numeric columns.
  if(ncol(sample_mat_numeric)!=ncol(sample_mat)){
    dplyr::bind_cols(equalized_sample_dt_complete,sample_mat_nonnumeric)->equalized_sample_dt_complete
  }
  
  return(equalized_sample_dt_complete)
}
