#' @title parallel_compare
#' @description This function performs a parallel comparison of the treatment effect of a control and treatment group.
#' @param long_ctrl_trt_df A long-format data frame with control and treatment ratio data.
#' @param interactive A logical value indicating whether to plot the results interactively (default is FALSE).
#' @param thresh the hurdle value, diving the binomial (zero or not) and the zero-inflated negative binomial (count) models. Values below this are not used in count models.
#' @param cap The outlier threshold.
#' @param debug A logical value indicating whether to enter the browser at the start of the function.
#' @param plot set to plot a heatmap of the results
#' @param cores The number of CPU cores to use for parallel processing (default is 1).
#' @return A data frame with the treatment effect comparison results.
#' @importFrom magrittr `%>%`
#' @importFrom ggplot2 geom_rect aes ggplot scale_fill_viridis_c scale_fill_viridis_c labs coord_flip facet_wrap theme
#' @examples
#' dontrun{
#' parallel_compare(long_ctrl_trt_df, mc.cores = 2, interactive = TRUE) 
#' }
#' @export
zinb_parallel_compare_glm=function(long_ctrl_trt_df,interactive=F,thresh=0,cap=Inf,debug=F,plot=F,cores=1,return_models=F){
  #browser()
  cell_ids=long_ctrl_trt_df %>% dplyr::filter(treatment=="trtlr") %>% dplyr::select(trtlr_cell_id) %>% unique() %>% dplyr::pull()
  if(debug){browser()}
  if(cores==1){
    zi_mods<-pbapply::pblapply(1:length(cell_ids),FUN=function(i) {
      glm(value==0 ~   treatment  ,data=long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt")  & trtlr_cell_id==cell_ids[i]) %>% dplyr::mutate(value=as.integer(value)) ,family=binomial)->zi_mod
      return(zi_mod)
    })
    zi_mods_tidy<-pbapply::pblapply(1:length(cell_ids),FUN=function(i) {
      broom::tidy(zi_mods[[i]]) %>% dplyr::mutate(model="zero",cell_id=cell_ids[i])-> zi_mod_tidy
      zi_mod_tidy$converged=zi_mods[[i]]$converged
      return(zi_mod_tidy)
    }) %>% do.call(rbind,.)
    nb_mods<-pbapply::pblapply(1:length(cell_ids),FUN=function(i) {
      tryCatch({
      MASS::glm.nb(value ~   treatment  ,data=long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt")  & trtlr_cell_id==cell_ids[i] & value>thresh) %>% dplyr::mutate(value=as.integer(value)) )->nb_mod
        return(nb_mod)
},error=function(e){warning(paste0(cell_ids[i]," failed to converge"));return();})
      
    })
    if(return_models){return(list(nb_mods=nb_mods,zi_mods=zi_mods))}
    #browser()
    nb_mods_tidy<-pbapply::pblapply(1:length(cell_ids),FUN=function(i) {
      if(is.null(nb_mods[[i]])){return(
        data.frame(term=c("(Intercept)","treatmenttrt"),estimate=c(0,0),std.error=c(1,1),statistic=c(0,0),p.value=c(1,1),converged=FALSE,model="trunc_nb",cell_id=cell_ids[i])
      )}
      tryCatch({
      broom::tidy(nb_mods[[i]])  %>% dplyr::mutate(model="trunc_nb",cell_id=cell_ids[i])-> nb_mod_tidy
        nb_mod_tidy$converged=nb_mods[[i]]$converged
      return(nb_mod_tidy)
      },error=function(e){warning(paste0(cell_ids[i]," failed to converge"));return(data.frame());})
    }) %>% do.call(rbind,.) #%>% dplyr::group_by(cell_id) %>% dplyr::mutate(cell_area_frac=cell_area/sum(cell_area)) %>% dplyr::ungroup() #,mc.cores=min(length(cell_ids),parallel::detectCores()) #
    dplyr::bind_rows(zi_mods_tidy,nb_mods_tidy) %>% dplyr::mutate(exp_est=exp(estimate),cell_area=get_cell_area(cell_id))   -> complete_model_results
    
    #lasso_results %>% dplyr::arrange(-exp_est)
    complete_model_results %>% dplyr::select(cell_id,cell_area) %>% dplyr::distinct() %>% dplyr::pull(cell_area) %>% sum()->total_cell_area
    complete_model_results %>% dplyr::mutate(cell_area_frac=cell_area/total_cell_area) %>% dplyr::select(cell_id,cell_area,cell_area_frac,dplyr::everything()) %>% dplyr::arrange(-exp_est) -> complete_model_results_with_cell_area
  }
  if(cores>1){
    
    zi_mods<-pbmcapply::pbmclapply(1:length(cell_ids),FUN=function(i) {
      glm(value==0 ~   treatment  ,data=long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt")  & trtlr_cell_id==cell_ids[i]) %>% dplyr::mutate(value=as.integer(value)) ,family=binomial)->zi_mod
      return(zi_mod)
    },mc.cores=cores)
    zi_mods_tidy<-pbmcapply::pbmclapply(1:length(cell_ids),FUN=function(i) {
      broom::tidy(zi_mods[[i]]) %>% dplyr::mutate(model="zero",cell_id=cell_ids[i])-> zi_mod_tidy
      zi_mod_tidy$converged=zi_mods[[i]]$converged
      return(zi_mod_tidy)
    },mc.cores=cores) %>% do.call(rbind,.)
    
    nb_mods<-pbmcapply::pbmclapply(1:length(cell_ids),FUN=function(i) {
      tryCatch({
        MASS::glm.nb(value ~   treatment  ,data=long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt")  & trtlr_cell_id==cell_ids[i] & value>thresh) %>% dplyr::mutate(value=as.integer(value)) )->nb_mod
        return(nb_mod)
      },error=function(e){warning(paste0(cell_ids[i]," failed to converge"));return();})
      
    },mc.cores=cores)
    if(return_models){return(list(nb_mods=nb_mods,zi_mods=zi_mods))}
    nb_pred_df_complete<-pbmcapply::pbmclapply(1:length(cell_ids),FUN=function(i) {
      tryCatch({
        long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt")  & trtlr_cell_id==cell_ids[i]) %>% dplyr::mutate(value=as.integer(value)) ->nb_pred_df
        predict(nb_mods[[i]],newdata=nb_data_df,type="response")  -> nb_pred_df$pred
        return(nb_pred_df)
      },error=function(e){warning(paste0(cell_ids[i]," failed to converge"));return(data.frame());})
      
    },mc.cores=cores) %>% do.call(rbind,.)
    
    nb_pred_df_ctrl<-pbmcapply::pbmclapply(1:length(cell_ids),FUN=function(i) {
      tryCatch({
        long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl")  & trtlr_cell_id==cell_ids[i]) %>% dplyr::mutate(value=as.integer(value)) ->nb_pred_df
        predict(nb_mods[[i]],newdata=nb_data_df,type="response")  -> nb_pred_df$pred
        return(nb_pred_df)
      },error=function(e){warning(paste0(cell_ids[i]," failed to converge"));return(data.frame());})
      
    },mc.cores=cores) %>% do.call(rbind,.)
    nb_pred_df_trt<-pbmcapply::pbmclapply(1:length(cell_ids),FUN=function(i) {
      tryCatch({
        long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("trt")  & trtlr_cell_id==cell_ids[i]) %>% dplyr::mutate(value=as.integer(value)) ->nb_pred_df
        predict(nb_mods[[i]],newdata=nb_data_df,type="response")  -> nb_pred_df$pred
        return(nb_pred_df)
      },error=function(e){warning(paste0(cell_ids[i]," failed to converge"));return(data.frame());})
      
    },mc.cores=cores) %>% do.call(rbind,.)
    nb_pred_df_complete=rbind(nb_pred_df_ctrl,nb_pred_df_trt)
    ggplot(data=nb_pred_df_complete,aes(x=start1,y=start2,color=value)) + geom_rect() + scale_color_viridis_c() + facet_wrap(~treatment,scales="free")
    #browser()
    nb_mods_tidy<-pbmcapply::pbmclapply(1:length(cell_ids),FUN=function(i) {
      if(is.null(nb_mods[[i]])){return(
        data.frame(term=c("(Intercept)","treatmenttrt"),estimate=c(0,0),std.error=c(1,1),statistic=c(0,0),p.value=c(1,1),converged=FALSE,model="trunc_nb",cell_id=cell_ids[i])
        )}
      tryCatch({
        broom::tidy(nb_mods[[i]])  %>% dplyr::mutate(model="trunc_nb",cell_id=cell_ids[i])-> nb_mod_tidy
        nb_mod_tidy$converged=nb_mods[[i]]$converged
        return(nb_mod_tidy)
      },error=function(e){warning(paste0(cell_ids[i]," failed to converge"));
        return(
          data.frame(term=c("(Intercept)","treatmenttrt"),estimate=c(0,0),std.error=c(1,1),statistic=c(0,0),p.value=c(1,1),converged=FALSE,model="trunc_nb",cell_id=cell_ids[i])
          )
        ;})   },mc.cores=cores) %>% do.call(rbind,.) #%>% dplyr::group_by(cell_id) %>% dplyr::mutate(cell_area_frac=cell_area/sum(cell_area)) %>% dplyr::ungroup() #,mc.cores=min(length(cell_ids),parallel::detectCores()) #
    dplyr::bind_rows(zi_mods_tidy,nb_mods_tidy) %>% dplyr::mutate(exp_est=exp(estimate),cell_area=get_cell_area(cell_id))   -> complete_model_results    
    #lasso_results %>% dplyr::arrange(-exp_est)
    complete_model_results %>% dplyr::select(cell_id,cell_area) %>% dplyr::distinct() %>% dplyr::pull(cell_area) %>% sum()->total_cell_area
    complete_model_results %>% dplyr::mutate(cell_area_frac=cell_area/total_cell_area) %>% dplyr::select(cell_id,cell_area,cell_area_frac,dplyr::everything()) %>% dplyr::arrange(-exp_est) -> complete_model_results_with_cell_area

    }
  #fdr adjust:
  #browser()
  # (
  #   complete_model_results_with_cell_area %>% dplyr::filter(start1>start2)  %>% dplyr::mutate(plotted_value=CNVScope::signedRescale((trt_wt_ratio))) 
  # ) %>% dplyr::bind_rows(
  #   
  #   complete_model_results_with_cell_area %>% dplyr::filter(start1<=start2)  %>% dplyr::mutate(plotted_value=ifelse(p.value!=0,-log(p.value),-log(min(complete_merged_glm_results$p.value[complete_merged_glm_results$p.value>0]))) %>% CNVScope::signedRescale())
  # ) %>% ggplot(aes(x=start1,y=start2,color=plotted_value,fill=plotted_value)) + geom_raster() + scale_color_viridis_c() + scale_fill_viridis_c() + theme_minimal() + labs(title="GLM results for trtlr cell id",x="Start1",y="Start2") -> pval_ratio_plot
  output_summary=complete_model_results_with_cell_area %>% tidyr::separate(cell_id,sep="_| ",into=c("chr2","start2","end2","chr1","start1","end1")) %>% dplyr::filter(model=="trunc_nb",term=="treatmenttrt") %>% dplyr::mutate(p.value=Rmpfr::pnorm(statistic,lower.tail = F)) %>% dplyr::mutate(padj=p.adjust(p.value,method="fdr"))
  #browser()
if(plot){
  print(output_summary  %>% ggplot(aes(xmin=as.numeric(start1),xmax=as.numeric(end1),ymin=as.numeric(start2),ymax=as.numeric(end2),fill=ifelse(padj<0.05,(abs(statistic)),0))) + geom_rect() + scale_fill_viridis_c(direction = 1,option="magma",oob=scales::squish)+ theme_minimal() + labs(title="truncated negative binomial GLMM Wald statistics (W=β/SE), by cell id",x="Start1",y="Start2")  + coord_flip() + theme(legend.position="bottom"))
  
  print(output_summary %>% dplyr::mutate(rate_ratio=as.numeric((exp(estimate)) ) ) %>% dplyr::arrange((rate_ratio))  %>% ggplot(aes(xmin=as.numeric(start1),xmax=as.numeric(end1),ymin=as.numeric(start2),ymax=as.numeric(end2),fill=rate_ratio)) + geom_rect() + scale_fill_viridis_c(option="magma")+ theme_minimal() + labs(title="truncated negative binomial GLMM rate ratio values, by cell id",x="Start1",y="Start2") + theme(legend.position="bottom"))}
  output_summary=complete_model_results_with_cell_area %>% tidyr::separate(cell_id,sep="_| ",into=c("chr2","start2","end2","chr1","start1","end1")) %>% dplyr::filter(term=="treatmenttrt") %>% dplyr::mutate(p.value=Rmpfr::pnorm(statistic,lower.tail = F)) %>% dplyr::mutate(padj=p.adjust(p.value,method="fdr")) %>% dplyr::arrange(term,-abs(statistic))  
  
  
  
  #if(interactive){plotly::ggplotly(pval_ratio_plot) %>% plotly::toWebGL()} else {print(pval_ratio_plot)}
  
  return(output_summary)
}
get_cell_area=function(cell_id) {
  (CNVScope::underscored_pos_to_GRanges(strsplit(cell_id," ")[[1]][[1]])@ranges@width)*(CNVScope::underscored_pos_to_GRanges(strsplit(cell_id," ")[[1]][[2]])@ranges@width) ->area
  return(area)
}