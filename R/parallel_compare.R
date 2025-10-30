#' @title parallel_compare
#' @description This function performs a parallel comparison of the treatment effect of a control and treatment group.
#' @param long_ctrl_trt_df A long-format data frame with control and treatment ratio data.
#' @param mc.cores The number of CPU cores to use for parallel processing (default is 1).
#' @param interactive A logical value indicating whether to plot the results interactively (default is FALSE).
#' @return A data frame with the treatment effect comparison results.
#' @importFrom magrittr `%>%`
#' @examples
#' dontrun{
#' parallel_compare(long_ctrl_trt_df, mc.cores = 2, interactive = TRUE) 
#' }
#' @export
parallel_compare=function(long_ctrl_trt_df,mc.cores=1,interactive=F){
  cell_ids=long_ctrl_trt_df %>% dplyr::filter(treatment=="trtlr") %>% dplyr::select(trtlr_cell_id) %>% unique() %>% dplyr::pull()
  if(mc.cores==1){
    zi_mods<-pbapply::pblapply(1:length(cell_ids),FUN=function(i) {
      glmmTMB::glmmTMB(value==0 ~   condition  ,data=long_ctrl_trt_df %>% dplyr::filter(condition %in% c("ctrl","trt") & treatment=="lasso" & trtlr_cell_id==cell_ids[i]) %>% dplyr::mutate(value=as.integer(value)) ,family=binomial)->zi_mod
      return(zi_mod)
    })
    zi_mods_tidy<-pbapply::pblapply(1:length(cell_ids),FUN=function(i) {
      broom.mixed::tidy(zi_mods[[i]]) %>% dplyr::mutate(model="zero",cell_id=cell_ids[i])-> zi_mod_tidy
      return(zi_mod_tidy)
    }) %>% do.call(rbind,.)
    nb_mods<-pbapply::pblapply(1:length(cell_ids),FUN=function(i) {
      glmmTMB::glmmTMB(value ~   condition  ,data=long_ctrl_trt_df %>% dplyr::filter(condition %in% c("ctrl","trt") & treatment=="lasso" & trtlr_cell_id==cell_ids[i] & value>thresh) %>% dplyr::mutate(value=as.integer(value)) ,family=glmmTMB::truncated_nbinom2())->nb_mod
      return(nb_mod)
    })
    nb_mods_tidy<-pbapply::pblapply(1:length(cell_ids),FUN=function(i) {
      broom.mixed::tidy(nb_mods[[i]])  %>% dplyr::mutate(model="trunc_nb",cell_id=cell_ids[i])-> nb_mod_tidy
      return(nb_mod_tidy)
    }) %>% do.call(rbind,.) #%>% dplyr::group_by(cell_id) %>% dplyr::mutate(cell_area_frac=cell_area/sum(cell_area)) %>% dplyr::ungroup() #,mc.cores=min(length(cell_ids),parallel::detectCores()) #
    dplyr::bind_rows(zi_mods_tidy,nb_mods_tidy) %>% dplyr::mutate(exp_est=exp(estimate),cell_area=get_cell_area(cell_id))   -> complete_model_results
    
    #lasso_results %>% dplyr::arrange(-exp_est)
    complete_model_results %>% dplyr::select(cell_id,cell_area) %>% dplyr::distinct() %>% dplyr::pull(cell_area) %>% sum()->total_cell_area
    complete_model_results %>% dplyr::mutate(cell_area_frac=cell_area/total_cell_area) %>% dplyr::select(cell_id,cell_area,cell_area_frac,dplyr::everything()) %>% dplyr::arrange(-exp_est) -> complete_model_results_with_cell_area
    lapply(1:length(cell_ids),function(i){
      cell_df=long_ctrl_trt_df %>% dplyr::filter(condition %in% c("ctrl","trt") & treatment=="lasso" & trtlr_cell_id==cell_ids[i] & value!=0) %>% dplyr::mutate(value=as.integer(value))
      predict(nb_mods[[i]],newdata=cell_df,type="response") %>% dplyr::mutate(cell_id=cell_ids[i]) %>% dplyr::mutate(model="trunc_nb") -> nb_mod_preds
      return(nb_mod_preds)
    }) %>% do.call(rbind,.)  %>% dplyr::mutate(cell_area=get_cell_area(cell_id)) %>% dplyr::mutate(cell_area_frac=cell_area/total_cell_area) %>% dplyr::select(cell_id,cell_area,cell_area_frac,model,cell_id,exp_est) %>% dplyr::arrange(-exp_est) -> nb_mod_preds_df
  }    
  if(mc.cores>1){
    
    zi_mods<-pbmcapply::pbmclapply(1:length(cell_ids),FUN=function(i) {
      glmmTMB::glmmTMB(value==0 ~   condition  ,data=long_ctrl_trt_df %>% dplyr::filter(condition %in% c("ctrl","trt") & treatment=="lasso" & trtlr_cell_id==cell_ids[i]) %>% dplyr::mutate(value=as.integer(value)) ,family=binomial)->zi_mod
      return(zi_mod)
    })
    zi_mods_tidy<-pbmcapply::pbmclapply(1:length(cell_ids),FUN=function(i) {
      broom.mixed::tidy(zi_mods[[i]]) %>% dplyr::mutate(model="zero",cell_id=cell_ids[i])-> zi_mod_tidy
      return(zi_mod_tidy)
    }) %>% do.call(rbind,.)
    nb_mods<-pbmcapply::pbmclapply(1:length(cell_ids),FUN=function(i) {
      glmmTMB::glmmTMB(value ~   condition  ,data=long_ctrl_trt_df %>% dplyr::filter(condition %in% c("ctrl","trt") & treatment=="lasso" & trtlr_cell_id==cell_ids[i] & value>thresh) %>% dplyr::mutate(value=as.integer(value)) ,family=glmmTMB::truncated_nbinom2())->nb_mod
      return(nb_mod)
    })
    nb_mods_tidy<-pbmcapply::pbmclapply(1:length(cell_ids),FUN=function(i) {
      broom.mixed::tidy(nb_mods[[i]])  %>% dplyr::mutate(model="trunc_nb",cell_id=cell_ids[i])-> nb_mod_tidy
      return(nb_mod_tidy)
    }) %>% do.call(rbind,.) #%>% dplyr::group_by(cell_id) %>% dplyr::mutate(cell_area_frac=cell_area/sum(cell_area)) %>% dplyr::ungroup() #,mc.cores=min(length(cell_ids),parallel::detectCores()) #
    dplyr::bind_rows(zi_mods_tidy,nb_mods_tidy) %>% dplyr::mutate(exp_est=exp(estimate),cell_area=get_cell_area(cell_id))   -> complete_model_results
     
    #lasso_results %>% dplyr::arrange(-exp_est)
    complete_model_results %>% dplyr::select(cell_id,cell_area) %>% dplyr::distinct() %>% dplyr::pull(cell_area) %>% sum()->total_cell_area
    complete_model_results %>% dplyr::mutate(cell_area_frac=cell_area/total_cell_area) %>% dplyr::select(cell_id,cell_area,cell_area_frac,dplyr::everything()) %>% dplyr::arrange(-exp_est) -> complete_model_results_with_cell_area
    # lapply(1:length(cell_ids),function(i){
    #   cell_df=long_ctrl_trt_df %>% dplyr::filter(condition %in% c("ctrl","trt") & treatment=="lasso" & trtlr_cell_id==cell_ids[i] & value!=0) %>% dplyr::mutate(value=as.integer(value))
    #   predict(nb_mods[[i]],newdata=cell_df,type="response") %>% dplyr::mutate(cell_id=cell_ids[i]) %>% dplyr::mutate(model="trunc_nb") -> nb_mod_preds
    #   return(nb_mod_preds)
    # }) %>% do.call(rbind,.)  %>% dplyr::mutate(cell_area=get_cell_area(cell_id)) %>% dplyr::mutate(cell_area_frac=cell_area/total_cell_area) %>% dplyr::select(cell_id,cell_area,cell_area_frac,model,cell_id,exp_est) %>% dplyr::arrange(-exp_est) -> nb_mod_preds_df    
    
  }
  #fdr adjust:
  #browser()
  (
    complete_model_results_with_cell_area %>% dplyr::filter(start1>start2)  %>% dplyr::mutate(plotted_value=CNVScope::signedRescale((trt_wt_ratio))) 
  ) %>% dplyr::bind_rows(
    
    complete_model_results_with_cell_area %>% dplyr::filter(start1<=start2)  %>% dplyr::mutate(plotted_value=ifelse(p.value!=0,-log(p.value),-log(min(complete_merged_glm_results$p.value[complete_merged_glm_results$p.value>0]))) %>% CNVScope::signedRescale())
  ) %>% ggplot(aes(x=start1,y=start2,color=plotted_value,fill=plotted_value)) + geom_raster() + scale_color_viridis_c() + scale_fill_viridis_c() + theme_minimal() + labs(title="GLM results for trtlr cell id",x="Start1",y="Start2") -> pval_ratio_plot
  if(interactive){plotly::ggplotly(pval_ratio_plot) %>% plotly::toWebGL()} else {print(pval_ratio_plot)}
  
  return(complete_merged_glm_results %>% dplyr::select(trtlr_cell_id,trt_wt_ratio,estimate,std.error,statistic,p.value) %>% dplyr::distinct() %>% dplyr::mutate(padj=p.adjust(Rmpfr::pnorm(estimate),method="fdr")))
}