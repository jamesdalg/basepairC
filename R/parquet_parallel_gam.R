#' Parallelized GAM Model Fitting for Parquet Datasets
#'
#' This function fits Generalized Additive Models (GAMs) to a parquet dataset, grouped by a specific cell identifier.
#' The computation is parallelized for efficiency, with options for multi-core processing and clearing previous results.
#'
#' @param parquet_file A file path to the parquet dataset to be used for the model fitting.
#' @param temp_dir A directory path to store temporary results. Defaults to `./tidy_objects`.
#' @param single_model_cores An integer specifying the number of CPU cores to use for each GAM model. Defaults to 1.
#' @param dt_cores An integer specifying the number of threads for `data.table` operations. Defaults to 1.
#' @param mc.cores An integer specifying the number of cores to use for parallel processing. If `mc.cores` is 1, `pbapply::pblapply()` is used; otherwise, `parallel::mclapply()` is used.
#' @param clear_results A logical value. If TRUE, clears any existing results in the `temp_dir` before processing. Defaults to FALSE.
#' @param cutpoint A value to filter the cell data frame, removing points below this mark.
#' @return A data.table containing the summary of fitted models for each cell, including coefficients, standard errors, AIC values, and convergence status.
#'
#' @details
#' The function processes parquet datasets in parallel by fitting a GAM for each unique `lasso_cell_id`. 
#' It applies corrections to each sample in the dataset using the `mnase_correct_by_sample` function before modeling the effect of treatment. 
#' The model is fitted using a Negative Binomial family with an identity link. 
#' Results for each `lasso_cell_id` are stored in the `temp_dir` as CSV files, and previously stored results can be cleared if `clear_results` is set to TRUE.
#'
#' Parallelization is handled by `pbapply::pblapply()` or `parallel::mclapply()` depending on the number of cores specified.
#'
#' @import mgcv
#' @import pbapply
#' @import parallel
#' @import data.table
#' @import dplyr
#' @import broom
#' @import arrow
#' @importFrom stats AIC
#'
#' @examples
#' \dontrun{
#'   results <- parquet_parallel_gam(parquet_file = "my_data.parquet", temp_dir = "./results", single_model_cores = 2, dt_cores = 4, mc.cores = 4)
#' }
#'
#' @export
parquet_parallel_gam=function(parquet_file,temp_dir="./tidy_objects",single_model_cores=1,dt_cores=1,mc.cores=1,clear_results=F,spec_cell_ids=NULL,cutpoint=0){
  pbapply_funs=list(pblapply=pbapply::pblapply, pbmcapply=pbmcapply::pbmclapply,mcl=parallel::mclapply)
  if(mc.cores==1){
    app_fun=pbapply_funs$pblapply
  }    else {
    app_fun=pbapply_funs$mcl
  }
  arrow::open_dataset(parquet_file)->parquet_cell_dataset
  dir.create(temp_dir,showWarnings = F,recursive = T)
  list.files(temp_dir)->tidy_obs
  if(clear_results){
    sapply(tidy_obs,function(x) file.remove(paste0(temp_dir,x)))
  }
  #parquet_cell_dataset %>% dplyr::filter(value>0) %>% dplyr::group_by(lasso_cell_id) %>%   dplyr::summarize(mean_value=mean(value)) %>% dplyr::arrange(-mean_value) %>% dplyr::collect() ->parquet_cell_id_mean_summary
  if(!is.null(spec_cell_ids)){
    cell_ids=spec_cell_ids
  } else{
  parquet_cell_dataset %>% dplyr::select(lasso_cell_id) %>% dplyr::distinct() %>% dplyr::collect() %>% unlist() %>% unique()->cell_ids
  }
  parquet_dataset %>% dplyr::filter(value>0) %>% dplyr::rename(cell_id=lasso_cell_id) %>% dplyr::group_by(cell_id) %>% dplyr::summarize(n_treatments=n_distinct(treatment)) %>% dplyr::filter(n_treatments>=2) %>% dplyr::collect()->cells_with_data
  #parquet_cell_dataset  %>% dplyr::group_by(lasso_cell_id) %>%   dplyr::summarize(n_samples=dplyr::n_distinct(sample),n_treat=dplyr::n_distinct(treatment))   %>% dplyr::collect() ->parquet_cell_id_treatment_summary
  options(single_model_cores=single_model_cores)
  setDTthreads(dt_cores)
  options(mc.cores=mc.cores)
  app_fun(1:length(cell_ids), function(i) {
    message_parallel(paste0("modeling cell:",cell_ids[i]))
    tryCatch({
      if(file.exists(paste0(temp_dir, cell_ids[i], ".csv"))){
        return(data.table::fread(paste0(temp_dir, cell_ids[i], ".csv")))
      }
      message_parallel(paste0("retreiving cell ",cell_ids[i]," ",i,"/",length(cell_ids)))
      cell_df <- parquet_cell_dataset %>% dplyr::filter(lasso_cell_id == cell_ids[i]) %>% dplyr::collect()
      
      
      cell_df_filtered=cell_df[cell_df$value>cutpoint,]
      cell_df_filtered[, treatment := as.factor(ifelse(grepl("dTAG", sample), "trt", "ctrl"))]
      cell_df_filtered[, sample:=as.factor(sample)]
      message_parallel(paste0(cell_ids[i],"number of samples",length(unique(cell_df_filtered$sample))))
      message_parallel(paste0(cell_ids[i],"number of treatments",length(unique(cell_df_filtered$treatment))))
      if((length(unique(cell_df_filtered$treatment))<2)){ #|(length(unique(cell_df_filtered$sample))<2)
        message_parallel(paste0(cell_ids[i]," has less than two samples or treatments",i))
        return(data.frame(term=c("(Intercept)","treatmenttrt"),estimate=c(0,0),std.error=c(1,1),statistic=c(0,0),p.value=c(1,1),model="nb",cell_id=cell_ids[i],AIC=0,converged=FALSE)
        )
      }
      rm(cell_df)
      cell_df_filtered %>% dplyr::rename(cell_id=lasso_cell_id)->cell_df_filtered
      cell_df_filtered[,total_reads:=sum(value),by=sample]
      message_parallel(paste0("correcting by sample ",cell_ids[i],i))
      mnase_correct_by_sample(cell_df_filtered)->cell_df_filtered_corrected
      message_parallel(paste0("modeling treatment ",cell_ids[i],i))
      system.time({gam_model_trt <- mgcv::gam(snr_filtered_value_gam ~   treatment,                      data = cell_df_filtered_corrected, family = mgcv::nb(link="identity"),control = list(nthread=single_model_cores))})
      gam_model_trt %>% broom::tidy(parametric=T) %>% dplyr::mutate(model="nb",cell_id=cell_ids[i])-> nb_mod_tidy
      nb_mod_tidy$AIC=AIC(gam_model_trt)
      nb_mod_tidy$converged=T
      data.table::fwrite(nb_mod_tidy, paste0(temp_dir, cell_ids[i], ".csv"))
      return(nb_mod_tidy)
    },error=function(e){warning(paste0(cell_ids[i]," failed to converge"));
      message_parallel(paste0(cell_ids[i]," failed to converge",i));
      nb_mod_tidy=data.frame(term=c("(Intercept)","treatmenttrt"),estimate=c(0,0),std.error=c(1,1),statistic=c(0,0),p.value=c(1,1),model="nb",cell_id=cell_ids[i],AIC=0,converged=FALSE)
      data.table::fwrite(nb_mod_tidy, paste0(temp_dir, cell_ids[i], ".csv"))
      return(
        nb_mod_tidy
      );})
  }) %>% rbindlist()->nb_mod_tidy_full
  return(nb_mod_tidy_full)
}
