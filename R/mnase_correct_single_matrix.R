#' Correct MNase Bias in a Single Matrix
#'
#' @description
#' This function corrects MNase bias in a single data matrix using a GAM model fit to the data. 
#' The function can use pre-computed bias matrices or a provided model and includes options for 
#' parallelization and model saving.
#'
#' @param data_mat A matrix containing the raw data (e.g., counts) to be corrected.
#' @param sample_num An optional integer specifying the number of samples to use from the data matrix. 
#' If `NULL`, all data will be used for fitting the model.
#' @param bias_mat An optional matrix of MNase bias values. If `NULL`, a default bias matrix will be loaded 
#' from the \code{basepairC} package.
#' @param snr_filter A logical value indicating whether to apply an SNR filter. Defaults to \code{FALSE}.
#' @param single_model_cores An integer specifying the number of cores to use for parallelization 
#' in model fitting. Defaults to the number of detected cores.
#' @param dtthreads An integer specifying the number of threads to use for \code{data.table} operations. 
#' Defaults to the number of detected cores.
#' @param correction_type A string specifying the type of correction to apply. Options are \code{"ratio"} 
#' or \code{"residual"}. Defaults to \code{"ratio"}.
#' @param pre_computed_model An optional pre-computed GAM model. If provided, the model will be used 
#' instead of fitting a new one.
#' @param save_model_loc An optional string specifying the file path to save the fitted GAM model. 
#' If \code{NULL}, the model will not be saved.
#'
#' @details
#' The function melts the input data and bias matrices to long format and merges them for analysis. 
#' A GAM model is fit to the merged data using \code{mgcv::gam}, with optional bootstrapping to 
#' determine starting values for the model coefficients. The data can be corrected using either a 
#' ratio or residual-based approach, depending on the \code{correction_type} argument.
#'
#' @return A matrix with the same dimensions as \code{data_mat}, containing the corrected values. 
#' The returned matrix includes the GAM model as an attribute \code{"gam_model"}.
#'
#' @examples
#' \dontrun{
#' # Correct a data matrix using the default bias matrix and all available cores
#' corrected_matrix <- mnase_correct_single_matrix(data_mat = my_data_matrix)
#'
#' # Use a subset of the data for model fitting
#' corrected_matrix <- mnase_correct_single_matrix(data_mat = my_data_matrix, sample_num = 10000)
#'
#' # Use a pre-computed model and save the new model to a file
#' corrected_matrix <- mnase_correct_single_matrix(
#'   data_mat = my_data_matrix, 
#'   pre_computed_model = my_precomputed_model, 
#'   save_model_loc = "path/to/save_model.rds"
#' )
#' }
#'
#' @importFrom parallel detectCores
#' @importFrom data.table setDT setDTthreads as.data.table merge.data.table setnames
#' @importFrom reshape2 melt
#' @importFrom mgcv gam nb
#' @importFrom dplyr sample_n
#' @export

mnase_correct_single_matrix=function(data_mat,sample_num=NULL,bias_mat=NULL,snr_filter=F,single_model_cores=parallel::detectCores(),dtthreads=parallel::detectCores(),correction_type="ratio",pre_computed_model=NULL,save_model_loc=NULL){
  
  if(is.null(bias_mat)){bias_mat=readRDS(system.file(package="basepairC","data/mnase_bias_6bp.rds"))}
  data.table::setDTthreads(dtthreads)
  #melt bias matrix to long form from wide matrix form.
  mnase_bias_matrix_melt <- reshape2::melt(bias_mat, id.vars = NULL,  value.name = "mnase") %>% data.table::as.data.table()
  #melt data matrix to long form from wide matrix form.
  data_mat_melt <- reshape2::melt(data_mat, id.vars = NULL,  value.name = "counts") %>% data.table::as.data.table()
  
  # Perform the join and renaming using data.table
  merged_dt <- data.table::merge.data.table(mnase_bias_matrix_melt, data_mat_melt, by = c("Var1", "Var2")) 
  data.table::setnames(merged_dt, old = c("Var1", "Var2", "mnase", "counts"), new = c("loc1", "loc2", "mnase", "counts"))
  data.table::setDT(merged_dt)
  merged_dt[, min_mnase := min(mnase)]
  merged_dt[, mnase := ztrans(log(mnase + min_mnase))]
  if(is.null(pre_computed_model)){
    starting_values <- bootstrap_gam_parallel(
      data = merged_dt, 
      formula = counts ~ mnase, 
      family = mgcv::nb(),
      fraction = 0.01, 
      n_bootstraps = 5
    )
    #    mgcv::gam(counts ~ mnase,data=merged_dt %>% dplyr::sample_n(1e5),family=mgcv::nb(),control=list(nthread=parallel::detectCores(),trace=T,optim=list(method="BFGS")),start=starting_values$avg_coefs)->mnase_gam_subsample
    if(is.null(sample_num)){
      mgcv::gam(counts ~ mnase,data=merged_dt ,family=mgcv::nb(),nthreads = parallel::detectCores(),control=list(nthread=parallel::detectCores(),trace=T,optim=list(method="BFGS")),start=starting_values$avg_coefs)->mnase_gam
    } else{
      mgcv::gam(counts ~ mnase,data=merged_dt %>% dplyr::sample_n(sample_num),family=mgcv::nb(),nthreads = parallel::detectCores(),control=list(nthread=parallel::detectCores(),trace=T,optim=list(method="BFGS")),start=starting_values$avg_coefs)->mnase_gam
    }
  } else{ mnase_gam=pre_computed_model}
  if(!is.null(save_model_loc)){saveRDS(mnase_gam,save_model_loc) }
  #mgcv::bam(counts ~ mnase,data=merged_dt,family=mgcv::nb(),discrete=T,nthreads=single_model_cores,samfrac=0.1)->mnase_bam
  print(summary(mnase_gam))
  print("mnase*beta_mnase distribution:")
  print((coef(mnase_gam)["mnase"]*(merged_dt$mnase)  ) %>% summary())
  if(correction_type=="ratio"){
    merged_dt$corrected=merged_dt$counts/exp(coef(mnase_gam)["(Intercept)"] + coef(mnase_gam)["mnase"] * merged_dt$mnase )
  }
  if(correction_type=="residual"){
    merged_dt$corrected=merged_dt$counts-exp(coef(mnase_gam)["(Intercept)"] + coef(mnase_gam)["mnase"] * merged_dt$mnase) 
    merged_dt$corrected=ifelse(merged_dt$corrected<0,0,merged_dt$corrected)
  }
  matrix(merged_dt$corrected,nrow=nrow(data_mat),ncol=ncol(data_mat),byrow=F)->corrected_mat
  colnames(corrected_mat)=colnames(data_mat)
  rownames(corrected_mat)=rownames(data_mat)
  attr(corrected_mat,"gam_model")=mnase_gam
  return(corrected_mat)
}
