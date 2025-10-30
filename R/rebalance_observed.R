#' Rebalance Observed Samples
#'
#' This function equalizes sample coverage by processing numeric columns in the input sample matrix.
#' It separates numeric and non-numeric columns, calculates total counts and non-zero counts, and 
#' rebalances the observed samples based on the specified filter type ("sorted" or "ranked").
#'
#' @param sample_mat A matrix or data frame containing the samples to be rebalanced. Numeric columns will be processed, while non-numeric columns will be kept intact.
#' @param filter_type A character string indicating the type of filtering to be applied. 
#' Can be either "sorted" (default, meaning to sort by the proportion of observed samples at each location followed by counts)
#'  or "ranked" (which takes the product of the two aforementioned factors and ranks them). 
#' 
#' @return A data table containing the rebalanced samples, with numeric columns equalized according to the selected filter type.
#'
#' @details
#' The function first identifies numeric columns in the input matrix. It then separates numeric and non-numeric columns,
#' calculates the total counts and the total non-zero counts for the numeric columns, and determines the minimum non-zero count.
#' Based on the filter type, the function either sorts or ranks the data to filter out excess reads, ensuring that the resulting
#' data table has equalized sample coverage across all numeric columns.
#'
#' 
#' @examples
#' \dontrun{
#' sample_data <- data.table::data.table(A = c(10, 0, 5), B = c(3, 5, 0), C = c("Type1", "Type2", "Type1"))
#' rebalanced_data <- rebalance_observed(sample_data, filter_type = "sorted")
#' }
#'
#' @export
rebalance_observed=function(sample_mat,filter_type="sorted"){
  message_parallel("Equalizing sample coverage...")
  message_parallel("identifying numeric columns")
  sapply(data.table::as.data.table(sample_mat),is.numeric)->numeric_columns
  message_parallel("separating numeric and non-numeric columns")
  #browser()
  sample_mat=data.table::as.data.table(sample_mat)
  sample_mat[,!..numeric_columns]->sample_mat_nonnumeric
  sample_mat[,..numeric_columns]->sample_mat_numeric
  n_samples=ncol(sample_mat_numeric)
  message_parallel("calculating total counts")
  total_counts=Rfast::colsums(as.matrix(sample_mat_numeric),parallel=0)
  message_parallel("creating nonzero matrix")
  sample_mat_nz=sample_mat_numeric
  sample_mat_nz[sample_mat_nz!=0]<-1
  message_parallel("calculating total nonzero counts")
  total_nz_counts=Rfast::colsums(as.matrix(sample_mat_nz),parallel=0)
  min(total_nz_counts)->min_total_nz_counts
  max_cores=parallel::detectCores() 
  outer_cores=min(ncol(sample_mat_numeric),max_cores)
  inner_cores=max(max_cores-outer_cores*n_samples,1)
  parallel::mclapply(1:ncol(sample_mat_numeric),function(i){
    message_parallel(paste0(" Column: ",i))
    sample_col=as.matrix(sample_mat_numeric)[,i]
    current_condition=conditions[i]
    cond_indicies=which(conditions==current_condition)
    cond_mat=sample_mat_nz[,..cond_indicies]
    message_parallel(paste0(i,"calculating the proportion of samples that have observed each loci in the matrix"))
    cond_prop_observed=Rfast::rowsums(as.matrix(cond_mat))/ncol(cond_mat)
    sample_col=as.matrix(sample_mat_numeric)[,i]
    message_parallel(paste0(i,"calculating the probability of selecting a read based on total counts"))
    prob_selected_read_prop=sample_col/sum(sample_col)
    if(filter_type=="ranked"){
      rank(cond_prop_observed*prob_selected_read_prop)->current_col_ranked
      sample_df=data.table::data.table(sample_col=sample_col,rank=current_col_ranked,sample_mat_rownum=1:nrow(sample_mat_numeric),prop_obs=cond_prop_observed,prob_selected_read_prop=prob_selected_read_prop)
      sample_df %>% dplyr::arrange(-rank) %>% head(min_total_nz_counts)->sample_df_filtered}
    if(filter_type=="sorted"){
      message_parallel(paste0(i,"sorting data"))
      sample_df=data.table::data.table(sample_col=sample_col,sample_mat_rownum=1:nrow(sample_mat_numeric),prop_obs=cond_prop_observed,prob_selected_read_prop=prob_selected_read_prop)
      sample_df[sample_df$sample_col!=0,] %>%  dplyr::arrange(-prop_obs,prob_selected_read_prop) %>% head(min_total_nz_counts)->sample_df_filtered
    }
    message_parallel(paste0(i,"filtering data"))
    sample_col_filtered=sample_col
    #sample_col_filtered[-sample_df_filtered$sample_mat_rownum]<-0
    sample_col_filtered[setdiff(1:length(sample_col),sample_df_filtered$sample_mat_rownum)]<-0
    sample_col_filtered=data.table::data.table(sample_col_filtered)
    colnames(sample_col_filtered)=colnames(sample_mat_numeric)[i]
    message_parallel(paste0(i,"finished"))
    return(sample_col_filtered)
  },mc.cores=outer_cores) %>% do.call(cbind,.)   -> equalized_sample_mat #%>% Rfast::rowmeans()
  data.table::data.table(equalized_sample_mat)->equalized_sample_dt

  if(ncol(sample_mat_numeric)!=ncol(sample_mat)){
    dplyr::bind_cols(equalized_sample_dt,sample_mat_nonnumeric)->equalized_sample_dt_complete
  } else{
    equalized_sample_dt_complete=equalized_sample_dt
  }
  
  return(equalized_sample_dt_complete)
  
}

