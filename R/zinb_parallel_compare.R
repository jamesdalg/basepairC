#' @title parallel_compare
#' @description This function performs a parallel comparison of the treatment effect of a control and treatment group.
#' @param long_ctrl_trt_df A long-format data frame with control and treatment ratio data.
#' @param mc.cores The number of CPU cores to use for parallel processing (default is 1).
#' @param interactive A logical value indicating whether to plot the results interactively (default is FALSE).
#' @param thresh the hurdle value, diving the binomial (zero or not) and the zero-inflated negative binomial (count) models. Values below this are not used in count models.
#' @param cap The outlier threshold.
#' @param influence Should influence measures be calculated?
#' @param return_models Should the models be returned?
#' @param single_model_cores The number of cores to use for single models.
#' @return A data frame with the treatment effect comparison results.
#' @importFrom magrittr `%>%`
#' @examples
#' dontrun{
#' parallel_compare(long_ctrl_trt_df, mc.cores = 2, interactive = TRUE) 
#' }
#' @export
zinb_parallel_compare=function(long_ctrl_trt_df,mc.cores=1,interactive=F,thresh=0,cap=Inf,fam=glmmTMB::truncated_nbinom2(),influence=F,return_models=F,offset=T,skip_zi=T,total_reads_ctrl=NULL,total_reads_trt=NULL,
                               glm_fn="glmmTMB",bias=F,bootstrap=F,n_straps=1e3,control_args=glmmTMB::glmmTMBControl(optimizer=optim,collect=F,profile=T,optArgs=list(method="CG"),parallel=1),single_model_cores=1){
if(is.null(control_args)){control_args=glmmTMB::glmmTMBControl(optimizer=optim,collect=F,profile=T,optArgs=list(method="CG"),parallel=1)}
control_args$parallel=single_model_cores
#  single_model_cores=1,
    #browser()
  options(mc.cores=mc.cores)
  pbapply_funs=list(pblapply=pbapply::pblapply, pbmcapply=pbmcapply::pbmclapply)
  if(mc.cores==1){
    app_fun=pbapply_funs$pblapply
  }    else {
    app_fun=parallel::mclapply
  }
  if(!is.null(total_reads_ctrl)){
    long_ctrl_trt_df %>% dplyr::mutate(total_reads=ifelse(treatment=="ctrl",total_reads_ctrl,total_reads))->long_ctrl_trt_df
  }
  if(!is.null(total_reads_trt)){
    long_ctrl_trt_df %>% dplyr::mutate(total_reads=ifelse(treatment=="trt",total_reads_trt,total_reads))->long_ctrl_trt_df
  }
  
  cell_ids=long_ctrl_trt_df %>% dplyr::filter(treatment=="trtlr") %>% dplyr::select(trtlr_cell_id) %>% unique() %>% dplyr::pull()
  data.table::setDTthreads(threads=mc.cores)
  if(mc.cores>length(cell_ids)){mc.cores=length(cell_ids)}
  #if(mc.cores==1){
  #browser()
  #control_args=glmmTMB::glmmTMBControl(optimizer=optim,collect=F,profile=T,
                                       #optArgs=list(method="CG"),parallel=single_model_cores)# ,iter.max=1e3,eval.max=1e3 #,reltol=1e-1#,maxit=1e5
  browser()
  if(skip_zi==F){
    zi_mods<-app_fun(1:length(cell_ids),FUN=function(i) {
      long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt")  & trtlr_cell_id==cell_ids[i]) %>% dplyr::mutate(value=as.integer(value))->cell_df
      if(cell_df %>% dplyr::filter(value==0) %>% nrow() > 0){
        glmmTMB::glmmTMB(value==0 ~   treatment  ,data= cell_df,family=binomial,control = control_args
                         )->zi_mod
      } else {
        zi_mod=NULL
      }
      #glmmTMB::glmmTMB(value==0 ~   treatment  ,data= cell_df,family=binomial)->zi_mod
      message_parallel(paste0("finished zi mod",i," of ",length(cell_ids)))
      return(zi_mod)
    })
browser()
    zi_mods_tidy_l<-app_fun(1:length(cell_ids),FUN=function(i) {
      if(is.null(zi_mods[[i]])){return(
        data.frame(effect="fixed",component="cond",term=c("(Intercept)","treatmenttrt"),estimate=c(0,0),std.error=c(1,1),statistic=c(0,0),p.value=c(1,1),model="zero",cell_id=cell_ids[i],converged=FALSE,AIC=0)
      )}
      broom.mixed::tidy(zi_mods[[i]]) %>% dplyr::mutate(model="zero",cell_id=cell_ids[i])-> zi_mod_tidy
      zi_mod_tidy$AIC=AIC(zi_mods[[i]])
      zi_mod_tidy$converged=T
      return(zi_mod_tidy)
    }) #%>% data.table::rbindlist()
    #browser()
    zi_mods_tidy=do.call(dplyr::bind_rows,zi_mods_tidy_l)
  } else { 
    zi_mods_tidy=data.frame()
  }
  # cell_id cell_area cell_area_frac effect component         term   estimate  std.error  statistic       p.value    model converged          AIC   exp_est
  # 1  chr15_61983326_61983465 chr15_61983326_61983465    383161     0.01515152  fixed      cond  (Intercept) 3.19623003 0.05142569  62.152397  0.000000e+00 trunc_nb      TRUE   3120.67396 24.440217
  # 2  chr15_61985426_61985675 chr15_61985426_61985675    383161     0.01515152  fixed      cond  (Intercept) 2.94407468 0.04242188  69.399916  0.000000e+00 trunc_nb      TRUE  10051.16918 18.993080
  # 3  chr15_61987106_61987495 chr15_61987106_61987495    383161     0.01515152  fixed      cond  (Intercept) 2.80203403 0.08305314  33.737847 1.611598e-249 trunc_nb      TRUE  26953.29908 16.478130
  # 4  chr15_61983466_61983605 chr15_61983466_61983605    383161     0.01515152  fixed      cond  (Intercept) 2.77292469 0.04600343  60.276471  0.000000e+00 trunc_nb      TRUE   2733.21045 16.005376
  # 5  chr15_61983326_61983465 chr15_61983466_61983605    383161     0.01515152  fixed      cond  (Intercept) 2.58336025 0.04323763  59.747964  0.000000e+00 trunc_nb      TRUE   5931.68101 13.241558
  # 6  chr15_61986956_61987105 chr15_61986956_61987105    383161     0.01515152  fixed      cond  (Intercept) 2.44032274 0.06643104  36.734675 2.042501e-295 trunc_nb      TRUE   2149.53659 11.476744
  #i=which(cell_ids=="chr15_61984436_61985515 chr15_61985676_61987145")
  i=which(cell_ids=="chr15_61983326_61983465 chr15_61983326_61983465")
  browser()
    nb_mods<-app_fun(1:length(cell_ids),FUN=function(i) {
      #check if there are at least two obs per condition.
      tryCatch({
        cell_df=long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt")  & trtlr_cell_id==cell_ids[i] & value>thresh) %>% dplyr::mutate(value=as.integer(value))
        if(length(unique(cell_df$treatment))<2 ){
          message_parallel(paste0("less than two conditions for model",i," of ",length(cell_ids)))
          return();
        }
        if(min(table(droplevels(cell_df$treatment)))<2){
          message_parallel(paste0("at least one condition has less than two observations for model",i," of ",length(cell_ids)))
          return()
        }
      browser()
      glmm_formula="value ~   treatment"
      
      if(offset){  glmm_formula=paste0(glmm_formula," + offset(log(total_reads))")}
      #glmmTMB::glmmTMB(value ~   treatment + offset(log(total_reads)) + splines::bs(mnase_bias),data= cell_df,family=fam,control=control_args)->nb_mod
         
      #else{
        #glmmTMB::glmmTMB(value ~   treatment ,data= cell_df,family=fam,control=control_args)->nb_mod
        #}
      if(bias){
        glmm_formula=paste0(glmm_formula," + splines::bs(mnase_bias)")
      }
      
      glmmTMB::glmmTMB(as.formula(glmm_formula),data= cell_df,family=fam,control=control_args)->nb_mod
      if(bootstrap){
        nb_mod=bootstrap_glmm_nb_hurdle(cell_name=cell_ids[i],formula_text=glmm_formula,long_ctrl_trt_df,fam=glmmTMB::truncated_nbinom2(),n_straps=1e3,n.cores=parallel::detectCores()/4,strap_size=NULL,return_full_results=F)
      }
      message_parallel(paste0("convergence; finished model",i," of ",length(cell_ids)))
      return(nb_mod)
      },error=function(e){warning(paste0(cell_ids[i]," failed to converge",i));return();})
      message_parallel(paste0("no convergence; finished model",i," of ",length(cell_ids)))
    })
    if(return_models){return(list(nb_mods=nb_mods,zi_mods=zi_mods))}
    nb_mods_tidy<-app_fun(1:length(cell_ids),FUN=function(i) {
      if(is.null(nb_mods[[i]])){return(
        data.frame(effect="fixed",component="cond",term=c("(Intercept)","treatmenttrt"),estimate=c(0,0),std.error=c(1,1),statistic=c(0,0),p.value=c(1,1),model="trunc_nb",cell_id=cell_ids[i],converged=FALSE,AIC=0)
      )}
      if(class(nb_mods[[i]]) %in% "data.frame"){
        nb_mods[[i]] %>% dplyr::mutate(model="trunc_nb",cell_id=cell_ids[i],converged=T,AIC=0)->nb_mod_tidy
      }
      tryCatch({
      broom.mixed::tidy(nb_mods[[i]])  %>% dplyr::mutate(model="trunc_nb",cell_id=cell_ids[i],converged=TRUE)-> nb_mod_tidy
        nb_mod_tidy$AIC=AIC(nb_mods[[i]])
      return(nb_mod_tidy)
      },error=function(e){warning(paste0(cell_ids[i]," failed to converge"));
        message_parallel(paste0(cell_ids[i]," failed to converge",i));
        return(
        data.frame(effect="fixed",component="cond",term=c("(Intercept)","treatmenttrt"),estimate=c(0,0),std.error=c(1,1),statistic=c(0,0),p.value=c(1,1),model="trunc_nb",cell_id=cell_ids[i],converged=FALSE,AIC=0)
        );})
    }) %>% data.table::rbindlist() # %>% do.call(rbind,.) #%>% dplyr::group_by(cell_id) %>% dplyr::mutate(cell_area_frac=cell_area/sum(cell_area)) %>% dplyr::ungroup() #,mc.cores=min(length(cell_ids),parallel::detectCores()) #
    #browser()
    dplyr::bind_rows(zi_mods_tidy,nb_mods_tidy) %>% dplyr::mutate(exp_est=exp(estimate),cell_area=get_cell_area(cell_id))   -> complete_model_results
    
    #lasso_results %>% dplyr::arrange(-exp_est)
    complete_model_results %>% dplyr::select(cell_id,cell_area) %>% dplyr::distinct() %>% dplyr::pull(cell_area) %>% sum()->total_cell_area
    complete_model_results %>% dplyr::mutate(cell_area_frac=cell_area/total_cell_area) %>% dplyr::select(cell_id,cell_area,cell_area_frac,dplyr::everything()) %>% dplyr::arrange(-exp_est) -> complete_model_results_with_cell_area
  #}
  # if(mc.cores>1){
  #   
  #   zi_mods<-app_fun(1:length(cell_ids),FUN=function(i) {
  #     basepaiC::message_parallel(paste0("zi mod",i," of ",length(cell_ids)))
  #     glmmTMB::glmmTMB(value==0 ~   treatment  ,data=long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt")  & trtlr_cell_id==cell_ids[i]) %>% dplyr::mutate(value=as.integer(value)) ,family=binomial)->zi_mod
  #     basepaiC::message_parallel(paste0("zi mod",i," of ",length(cell_ids)))
  #     return(zi_mod)
  #   },mc.cores=mc.cores)
  #   zi_mods_tidy<-app_fun(1:length(cell_ids),FUN=function(i) {
  #     broom.mixed::tidy(zi_mods[[i]]) %>% dplyr::mutate(model="zero",cell_id=cell_ids[i])-> zi_mod_tidy
  #     zi_mod_tidy$AIC=AIC(zi_mod[[i]])
  #     return(zi_mods_tidy)
  #   },mc.cores=mc.cores) %>% data.table::rbindlist()
  #   #browser()
  #   nb_mods<-app_fun(1:length(cell_ids),FUN=function(i) {
  #     tryCatch({
  #     long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt")  & trtlr_cell_id==cell_ids[i] & value>thresh) %>% dplyr::mutate(value=as.integer(value)) ->long_ctrl_trt_df_cell
  #     glmmTMB::glmmTMB(value ~   treatment  ,data= long_ctrl_trt_df_cell,family=fam)->nb_mod
  #     return(nb_mod)
  #     },error=function(e){warning(paste0(cell_ids[i]," failed to converge"));return();})
  #   },mc.cores=mc.cores)
  #   if(return_models){return(list(nb_mods=nb_mods,zi_mods=zi_mods))}
  #   nb_mods_tidy<-app_fun(1:length(cell_ids),FUN=function(i) {
  #     if(is.null(nb_mods[[i]])){return(
  #       data.frame(term=c("(Intercept)","treatmenttrt"),component="cond",estimate=c(0,0),std.error=c(1,1),statistic=c(0,0),p.value=c(1,1),converged=FALSE,model="trunc_nb",cell_id=cell_ids[i])
  #     )}
  #     tryCatch({
  #     broom.mixed::tidy(nb_mods[[i]])  %>% dplyr::mutate(model="trunc_nb",cell_id=cell_ids[i])-> nb_mod_tidy
  #       nb_mod_tidy$AIC=AIC(nb_mods[[i]])
  #     return(nb_mod_tidy)
  #   },error=function(e){warning(paste0(cell_ids[i]," failed to converge"));return(data.frame());})
  #   },mc.cores=mc.cores) %>% data.table::rbindlist() #%>% dplyr::group_by(cell_id) %>% dplyr::mutate(cell_area_frac=cell_area/sum(cell_area)) %>% dplyr::ungroup() #,mc.cores=min(length(cell_ids),parallel::detectCores()) #
  #   #browser()
  #   dplyr::bind_rows(zi_mods_tidy,nb_mods_tidy) %>% dplyr::mutate(exp_est=exp(estimate),cell_area=get_cell_area(cell_id))   -> complete_model_results
  #   
  #   #lasso_results %>% dplyr::arrange(-exp_est)
  #   complete_model_results %>% dplyr::select(cell_id,cell_area) %>% dplyr::distinct() %>% dplyr::pull(cell_area) %>% sum()->total_cell_area
  #   complete_model_results %>% dplyr::mutate(cell_area_frac=cell_area/total_cell_area) %>% dplyr::select(cell_id,cell_area,cell_area_frac,dplyr::everything()) %>% dplyr::arrange(-exp_est) -> complete_model_results_with_cell_area
  # 
  #   }
  #fdr adjust:
  #browser()
  # (
  #   complete_model_results_with_cell_area %>% dplyr::filter(start1>start2)  %>% dplyr::mutate(plotted_value=CNVScope::signedRescale((trt_wt_ratio))) 
  # ) %>% dplyr::bind_rows(
  #   
  #   complete_model_results_with_cell_area %>% dplyr::filter(start1<=start2)  %>% dplyr::mutate(plotted_value=ifelse(p.value!=0,-log(p.value),-log(min(complete_merged_glm_results$p.value[complete_merged_glm_results$p.value>0]))) %>% CNVScope::signedRescale())
  # ) %>% ggplot(aes(x=start1,y=start2,color=plotted_value,fill=plotted_value)) + geom_raster() + scale_color_viridis_c() + scale_fill_viridis_c() + theme_minimal() + labs(title="GLM results for trtlr cell id",x="Start1",y="Start2") -> pval_ratio_plot
  #browser()
  output_summary=complete_model_results_with_cell_area %>% tidyr::separate(cell_id,sep="_| ",into=c("chr2","start2","end2","chr1","start1","end1")) %>% dplyr::mutate(p.value=Rmpfr::pnorm(abs(statistic),lower.tail = F)) %>% dplyr::mutate(padj=p.adjust(p.value,method="fdr")) %>% dplyr::filter(model=="trunc_nb",term=="treatmenttrt") 

  plot(output_summary  %>% ggplot(aes(xmin=as.numeric(start1),xmax=as.numeric(end1),ymin=as.numeric(start2),ymax=as.numeric(end2),fill=ifelse(padj<0.05,(abs(statistic)),0))) + geom_rect() + scale_fill_viridis_c(direction = 1,option="magma",oob=scales::squish)+ theme_minimal() + labs(title="truncated negative binomial GLMM Wald statistics (W=β/SE), by cell id",x="Start1",y="Start2")  ) + ggplot2::coord_equal() + ggplot2::theme(legend.position="bottom")
  
  plot(output_summary %>% dplyr::mutate(rate_ratio=as.numeric((exp(estimate)) ) ) %>% dplyr::arrange((rate_ratio))  %>% ggplot(aes(xmin=as.numeric(start1),xmax=as.numeric(end1),ymin=as.numeric(start2),ymax=as.numeric(end2),fill=rate_ratio)) + geom_rect() + scale_fill_viridis_c(option="magma")+ theme_minimal() + labs(title="truncated negative binomial GLMM rate ratio values, by cell id",x="Start1",y="Start2")) + ggplot2::coord_equal() + ggplot2::theme(legend.position="bottom")
  #if(interactive){plotly::ggplotly(pval_ratio_plot) %>% plotly::toWebGL()} else {print(pval_ratio_plot)}
  output_summary=complete_model_results_with_cell_area %>% tidyr::separate(cell_id,sep="_| ",into=c("chr2","start2","end2","chr1","start1","end1"),remove = F)  %>% dplyr::mutate(p.value=Rmpfr::pnorm(abs(statistic),lower.tail = F)) %>% dplyr::mutate(padj=p.adjust(p.value,method="fdr")) %>% dplyr::arrange(term,-abs(statistic))
  return(list(output_summary=output_summary,cell_ids=cell_ids,complete_model_results_with_cell_area=complete_model_results_with_cell_area))
}
get_cell_area=function(cell_id) {
  (CNVScope::underscored_pos_to_GRanges(strsplit(cell_id," ")[[1]][[1]])@ranges@width)*(CNVScope::underscored_pos_to_GRanges(strsplit(cell_id," ")[[1]][[2]])@ranges@width) ->area
  return(area)
}