#' Perform Bootstrapped GLMM with Negative Binomial Hurdle Model
#'
#' This function performs a bootstrapped Generalized Linear Mixed Model (GLMM) using a truncated negative binomial family. 
#' The function applies the model to a specified cell's data, fitting the model multiple times with resampled data 
#' to estimate the variability of the model's parameters.
#'
#' @param cell_name Character. The name of the cell for which the bootstrap analysis is to be performed. The dataset is filtered based on this value.
#' @param formula_text Character. A formula in text format specifying the model to be fitted. Default is `"value ~ treatment"`.
#' @param long_ctrl_trt_df DataFrame. A long-format data frame containing the control and treatment data.
#' @param fam Family. The family to be used in the GLMM model. Default is `glmmTMB::truncated_nbinom2()`.
#' @param n_straps Numeric. The number of bootstrap iterations. Default is `1000`.
#' @param n.cores Numeric. The number of cores to be used for parallel processing. Default is `parallel::detectCores()/4`.
#' @param strap_size Numeric. The size of the bootstrap samples. If `NULL`, the sample size is set to the number of rows in the filtered data. Default is `NULL`.
#' @param return_full_results Logical. If `TRUE`, returns both the full bootstrap results and the summarized results. If `FALSE`, only the summarized results are returned. Default is `FALSE`.
#'
#' @return A data frame containing the summarized bootstrap results, including estimates, standard errors, AIC values, and confidence intervals. If `return_full_results` is `TRUE`, a list is returned with both the full bootstrap results and the summarized results.
#'
#' @examples
#' \dontrun{
#'   result <- bootstrap_glmm_nb_hurdle(
#'     cell_name = "cell_1",
#'     formula_text = "value ~ treatment",
#'     long_ctrl_trt_df = my_data_frame,
#'     n_straps = 500,
#'     n.cores = 4,
#'     strap_size = 1000,
#'     return_full_results = TRUE
#'   )
#' }
#'
#' @export
bootstrap_glmm_nb_hurdle=function(cell_name,formula_text="value ~ treatment",long_ctrl_trt_df,fam=glmmTMB::truncated_nbinom2(),n_straps=1e3,n.cores=parallel::detectCores()/4,strap_size=NULL,return_full_results=F){
  #browser()
  #profvis::profvis({
  #start timing
  pbapply_funs=list(pblapply=pbapply::pblapply, pbmcapply=pbmcapply::pbmclapply,mcl=parallel::mclapply)
  if(n.cores==1){
    app_fun=pbapply_funs$pblapply
  }    else {
    app_fun=pbapply_funs$mcl
  }
  proc.time()->start_time
  
  long_ctrl_trt_df %>% dplyr::filter(trtlr_cell_id==cell_name) %>% 
    dplyr::mutate(value=as.integer(value)) %>% 
    dplyr::mutate(treatment=factor(treatment,levels=c("ctrl","trt")))->cell_df
  if(is.null(strap_size)){
    nrow(cell_df)->strap_size
  }
  message_parallel("Starting bootstrap for cell ",cell_name)
  #time to create the bootstrap samples
  subset_time=proc.time()-start_time
  #boot_straps =
  #  tibble::tibble(
  #    strap_number = 1:n_straps,
  #    strap_sample =pbmcapply::pbmclapply(strap_number,function(iteration) {
  #  set.seed(iteration);
  #  return(sample(x = 1:nrow(cell_df),size = 1000,replace=T))
  #},mc.cores=n.cores))
  #proc.time()-start_time->bootstrap_sample_time
  #  basepairC::message_parallel("time to create bootstrap samples:","\n","user:",bootstrap_sample_time["user.self"],"\n system:",bootstrap_sample_time["sys.self"],"\n",
  #                              "elapsed:",bootstrap_sample_time["elapsed"],"cell:",cell_name)
  
  #bootstrap_results = 
  
  options(mc.cores=n.cores)
  tidy_models = app_fun(1:n_straps, function(strap_num) {
    tryCatch({
      set.seed(strap_num);
      #sample(x = 1:nrow(cell_df),size = 1000,replace=T)
      basepairC::message_parallel(paste0("running strap",strap_num," of ", n_straps))
      glmmTMB::glmmTMB(as.formula(formula_text), data = cell_df[unlist(sample(x = 1:nrow(cell_df),size = strap_size,replace=T)),],family=fam)->mod
      basepairC::message_parallel(paste0("finished strap",strap_num," of ", n_straps))
      broom.mixed::tidy(mod)->tidy_model;
      basepairC::message_parallel(paste0("tidy strap",strap_num," of ", n_straps))
      tidy_model$aic=AIC(mod)
      basepairC::message_parallel(paste0("AIC strap",strap_num," of ", n_straps))
      tidy_model$cell_id=cell_name
      
      tidy_model$strap_number=strap_num
      return(tidy_model)
    },error=function(e) {return(data.table::data.table())})
  } )
  model_time=proc.time()-subset_time
  basepairC::message_parallel("time to fit models:","\n","user:",model_time["user.self"],"\n system:",model_time["sys.self"],"\n",
                              "elapsed:",model_time["elapsed"],"cell:",cell_name)
  
  #aic time
  #aic_time=proc.time()-tidy_time
  # basepairC::message_parallel("time to calculate AIC:","\n","user:",aic_time["user.self"],"\n system:",aic_time["sys.self"],"\n",
  #                             "elapsed:",aic_time["elapsed"],"cell:",cell_name)
  
  bootstrap_results = tidy_models %>% data.table::rbindlist() #%>% dplyr::mutate(aic=aic_values,cell_id=cell_name)
  #rbind_time
  rbind_time=proc.time()-model_time
  basepairC::message_parallel("time to rbind:","\n","user:",rbind_time["user.self"],"\n system:",rbind_time["sys.self"],"\n",
                              "elapsed:",rbind_time["elapsed"],"cell:",cell_name)
  #browser()
  tryCatch({   bootstrap_results %>% na.omit() %>%  
      dplyr::group_by(term,cell_id) %>% 
      dplyr::summarize(boot_est=mean(estimate),boot_se = sd(estimate),boot_stat=boot_est/boot_se,boot_p=Rmpfr::pnorm(abs(boot_stat),lower.tail = F),mean_aic=mean(aic),sd_aic=sd(aic),var_aic=var(aic)) %>% 
      dplyr::arrange(abs(boot_se)) %>% dplyr::mutate(ci_l=boot_est-boot_se,ci_r=boot_est+boot_se,family=as.character(fam$family)
      ) %>% as.data.frame() ->bootstrapped_summarized
  },error=function(e) {bootstrap_results=data.table::data.table();return(bootstrap_results)})
  #summarize time
  summarize_time=proc.time()-rbind_time
  basepairC::message_parallel("time to summarize:","\n","user:",summarize_time["user.self"],"\n system:",summarize_time["sys.self"],"\n",
                              "elapsed:",summarize_time["elapsed"],"cell:",cell_name)
  #return(bootstrap_results)
  #})
  if(return_full_results){
    return(list(bootstrap_results_full=bootstrap_results,bootstrap_results_summarized=bootstrapped_summarized))
  }
  #plot the estimate over time with a linear regression line, with a beta coefficient and p-value for beta of the regression of the bootstrapped estimates (to determine convergence).
  #browser()
  plot(ggplot2::ggplot(bootstrap_results %>% dplyr::mutate(termcomp=paste(term,component)),ggplot2::aes(x=strap_number,y=estimate)) + ggplot2::geom_point() + ggplot2::geom_smooth(method="lm") + ggplot2::labs(title=paste0("Bootstrapped estimates of treatment effect for cell ",cell_name),x="Bootstrap iteration",y="Estimate of treatment effect")    +
         ggpubr::stat_regline_equation(
           ggplot2::aes(label =  paste(..eq.label.., ..adj.rr.label.., paste0("Var(\beta)=",var(estimate)),sep = "~~~~")),
           formula = estimate ~ strap_number,
         )  + ggplot2::facet_wrap(~termcomp))
  return(bootstrapped_summarized)
  #})
}