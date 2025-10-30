#' Create a Parquet Dataset from a List of Matrices
#'
#' This function processes a list of matrices, applies preprocessing steps, merges with a melted bias data table, and writes the resulting data to a Parquet dataset. 
#' It supports parallel processing and low disk footprint options.
#' Only choose the low disk footprint option if you have limited disk space and have enough memory to store the matrices (typically 512GB would be enough to comfortably do this).
#'
#' @param mat_list A list of matrices to be processed, where each matrix is associated with a sample.
#' @param segments A segmentation object used for preprocessing the matrices.
#' @param bias_dt_melted A melted data.table containing bias matrices, to be merged with the processed matrices.
#' @param parquet_dir A directory path where the Parquet files will be written. The directory will be created if it does not exist.
#' @param final_parquet_file Path to the final Parquet dataset that will be created after processing all matrices.
#' @param low_disk_footprint Logical flag indicating whether to use a low disk footprint approach (default is \code{FALSE}).
#' @param ncores Number of cores to use for parallel processing (default is 1).
#'
#' @details
#' This function applies the following steps to each matrix in \code{mat_list}:
#' \enumerate{
#'   \item Preprocesses the matrix using the provided segmentation object.
#'   \item Converts the preprocessed matrix to a data.table.
#'   \item Adds sample metadata to the data.table.
#'   \item Removes rows with specific treatment conditions.
#'   \item Extracts chromosome, start, and end positions from locational information.
#'   \item Renames the treatment column to distinguish between control and treatment groups.
#'   \item Merges the processed matrix with the bias data.table based on genomic coordinates.
#'   \item Writes the merged data to a Parquet file in the specified directory.
#' }
#' The function writes individual Parquet files for each matrix and then combines them into a Parquet dataset.
#'
#' @return \code{NULL}. The processed matrices are written to Parquet files in the specified directory.
#'
#' @examples
#' \dontrun{
#' create_parquet_set(mat_list, segments, bias_dt_melted, "output/parquet_dir", "final_parquet_file", low_disk_footprint = TRUE, ncores = 4)
#' }
#' 
#' @importFrom data.table setDTthreads as.data.table merge.data.table setnames rbindlist
#' @importFrom pbapply pblapply
#' @importFrom arrow write_parquet open_dataset
#' @export
create_parquet_set=function(mat_list,segments,bias_dt_melted,parquet_dir,final_parquet_file,low_disk_footprint=F,ncores=1,mc.cores=1,verbose=F) {
  data.table::setDTthreads(ncores)
  options(mc.cores=mc.cores)
  pbapply_funs=list(pblapply=pbapply::pblapply, pbmcapply=pbmcapply::pbmclapply)
  if(mc.cores==1){
    app_fun=pbapply_funs$pblapply
  }    else {
    app_fun=parallel::mclapply
  }
  #browser()
  app_fun(1:length(mat_list_no_scn_ds), function(i) {
    already_written=(!low_disk_footprint)&(file.exists(paste0(parquet_dir,"sample_mat_list_ds_df_ctcf_",attr(mat_list[[i]],"sample_name"),".parquet")))
    if(already_written){
      return(data.table())
    }
    #process the matrix
    preprocess_mat(fn = mat_list[[i]],verbose=T,dsfactor=1,segments_orig_obj = segments,segments_lasso_obj = segments,n_cores = ncores)->sample_mat_list_ds
    if(verbose){message_parallel("Processed matrix ",i," of ",length(mat_list))}
    #convert to data.table
    sample_mat_list_ds %>% data.table::as.data.table()->sample_mat_list_ds
    if(verbose){message_parallel("converted to data.table ",i," of ",length(mat_list))}
    #set sample name
    sample_mat_list_ds$sample=attr(mat_list[[i]],"sample_name")
    if(verbose){message_parallel("assigned_sample_name ",i," of ",length(mat_list))}
    #remove the lasso lines
    sample_mat_list_ds[treatment=="orig"]->sample_mat_list_ds
    if(verbose){message_parallel("removed_lasso_lines ",i," of ",length(mat_list))}
    #adding chr,start,end columns from loc1 and loc2
    sample_mat_list_ds[, c("chr1", "start1", "end1") := tstrsplit(loc1, split = "_| ")]
    if(verbose){message_parallel("split loc1 ",i," of ",length(mat_list))}
    sample_mat_list_ds[, c("chr2", "start2", "end2") := tstrsplit(loc2, split = "_| ")]
    if(verbose){message_parallel("split loc2 ",i," of ",length(mat_list))}
    sample_mat_list_ds[, `:=`(start1 = as.numeric(start1), start2 = as.numeric(start2))]
    if(verbose){message_parallel("converted to numeric ",i," of ",length(mat_list))}
    #rename the treatment column to be the one with ctrl and trt.
    sample_mat_list_ds[,treatment := NULL]
    data.table::setnames(sample_mat_list_ds,"treatment_lasso","treatment")
    sample_mat_list_ds[,treatment:=ifelse(grepl("WT",attr(mat_list[[i]],"sample_name")),"ctrl","trt")]
    if(verbose){message_parallel("reset treatment ",i," of ",length(mat_list))}
    #now, merge on the bias matrices
    data.table::merge.data.table(sample_mat_list_ds,bias_dt_melted,by=c("chr1","start1","end1","chr2","start2","end2"),all.x=T)->sample_mat_list_ds_merged
    if(verbose){message_parallel("merged bias ",i," of ",length(mat_list))}
    #if parquet directory doesn't exist, create it.
    dir.create(parquet_dir,showWarnings = F,recursive = T)
    if(verbose){message_parallel("created parquet directory ",i," of ",length(mat_list))}
    #write to parquet dataset
    if(low_disk_footprint) {return(sample_mat_list_ds_merged)} else {
      if(verbose){message_parallel("writing parquet file ",i," of ",length(mat_list))}
    arrow::write_parquet(sample_mat_list_ds_merged,paste0(parquet_dir,"sample_mat_list_ds_df_ctcf_",attr(mat_list[[i]],"sample_name"),".parquet"))
      if(verbose){message_parallel("written parquet file! ",i," of ",length(mat_list))}
    return(data.table()) }
  }) %>% data.table::rbindlist()->sample_mat_list_ds_df_full
  if(low_disk_footprint) {
    arrow::write_parquet(sample_mat_list_ds_df_full,final_parquet_file)
  } else {
    sample_mat_list_ds_df_parquet <- open_dataset(parquet_dir, format = "parquet")
    if(verbose){message_parallel("writing final parquet file ", final_parquet_file)}
    arrow::write_parquet(sample_mat_list_ds_df_parquet,final_parquet_file)
  }

  
}