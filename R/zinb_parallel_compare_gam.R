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
#' @param link The link function for the negative binomial model.
#' @param offset A logical value indicating whether to include an offset term in the model.
#' @param skip_zi A logical value indicating whether to skip the zero-inflated model.
#' @param total_reads_ctrl The total number of reads in the control group.
#' @param total_reads_trt The total number of reads in the treatment group.
#' @param mnase A logical value indicating whether to include MNase bias correction.
#' @param bootstrap A logical value indicating whether to perform bootstrapping.
#' @param n_straps The number of bootstrap samples.
#' @param mnase_mat A numeric matrix representing the MNase bias matrix.
#' @param fam family object from the mgcv package.
#' @param ridge Apply ridge penalty to treatment? Default FALSE (takes a bit longer).
#' @param lasso1d_filter Apply 1D LASSO filter to covariates? Default TRUE.
#' @param restricted_mnase Apply MNase bias correction only to the cell of interest? Default TRUE.
#' @param base_formula The base formula for the model.
#' @param robust_cov Use robust covariance matrix for p-values? Default FALSE.
#' @param plot_direction Direction for plotting. 1 for treatment > control, -1 for control > treatment.
#'  Changing this value to false models the entire matrix for MNase, which can be slow and memory-intensive.
#' @return A data frame with the treatment effect comparison results.
#' @importFrom magrittr `%>%`
#' @examples
#' dontrun{
#' parallel_compare(long_ctrl_trt_df, mc.cores = 2, interactive = TRUE) 
#' }
#' @export
zinb_parallel_compare_gam=function(long_ctrl_trt_df,mc.cores=1,interactive=F,thresh=0,cap=Inf,fam="nb",influence=F,return_models=F,offset=T,skip_zi=T,total_reads_ctrl=NULL,total_reads_trt=NULL,
                               mnase=F,bootstrap=F,n_straps=1e3,mnase_mat=NULL,single_model_cores=1,link="log",extra_formula_terms="",ridge=F,lasso1d_filter=F,restricted_mnase=T,base_formula="value ~   treatment",robust_cov=F,plot_direction=1){
#browser()
if(lasso1d_filter){
  message_parallel("performing 1D lasso filter");
  lasso_1d_filter(long_ctrl_trt_df,cap = cap,thresh=thresh,mc.cores = mc.cores)->dropped_covariates}
  #browser()
#if(is.null(control_args)){control_args=gamTMB::gamTMBControl(optimizer=optim,collect=F,profile=T,optArgs=list(method="CG"),parallel=1)}
#control_args$parallel=single_model_cores
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
  #control_args=gamTMB::gamTMBControl(optimizer=optim,collect=F,profile=T,
                                       #optArgs=list(method="CG"),parallel=single_model_cores)# ,iter.max=1e3,eval.max=1e3 #,reltol=1e-1#,maxit=1e5
  #browser()
  if(skip_zi==F){
    zi_mods<-app_fun(1:length(cell_ids),FUN=function(i) {
      long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt")  & trtlr_cell_id==cell_ids[i]) %>% dplyr::mutate(value=as.integer(value))->cell_df
      if(cell_df %>% dplyr::filter(value==0) %>% nrow() > 0){
        gam(value==0 ~   treatment  ,data= cell_df,family=binomial
                         )->zi_mod
      } else {
        zi_mod=NULL
      }
      #gamTMB::gamTMB(value==0 ~   treatment  ,data= cell_df,family=binomial)->zi_mod
      message_parallel(paste0("finished zi mod",i," of ",length(cell_ids)))
      return(zi_mod)
    })
#browser()
    zi_mods_tidy_l<-app_fun(1:length(cell_ids),FUN=function(i) {
      if(is.null(zi_mods[[i]])){return(
        data.frame(effect="fixed",component="cond",term=c("(Intercept)","treatmenttrt"),estimate=c(0,0),std.error=c(1,1),statistic=c(0,0),p.value=c(1,1),model="zero",cell_id=cell_ids[i],converged=FALSE,AIC=0)
      )}
      broom.mixed::tidy(zi_mods[[i]],parametric=T) %>% dplyr::mutate(model="zero",cell_id=cell_ids[i])-> zi_mod_tidy
      
      zi_mod_tidy$AIC=AIC(zi_mods[[i]])
      zi_mod_tidy$converged=T
      return(zi_mod_tidy)
    }) #%>% data.table::rbindlist()
    #browser()
    zi_mods_tidy=do.call(dplyr::bind_rows,zi_mods_tidy_l)
  } else { 
    zi_mods_tidy=data.frame()
    zi_mods=as.list(rep(NA,length(cell_ids)))
    
  }
  # cell_id cell_area cell_area_frac effect component         term   estimate  std.error  statistic       p.value    model converged          AIC   exp_est
  # 1  chr15_61983326_61983465 chr15_61983326_61983465    383161     0.01515152  fixed      cond  (Intercept) 3.19623003 0.05142569  62.152397  0.000000e+00 nb      TRUE   3120.67396 24.440217
  # 2  chr15_61985426_61985675 chr15_61985426_61985675    383161     0.01515152  fixed      cond  (Intercept) 2.94407468 0.04242188  69.399916  0.000000e+00 nb      TRUE  10051.16918 18.993080
  # 3  chr15_61987106_61987495 chr15_61987106_61987495    383161     0.01515152  fixed      cond  (Intercept) 2.80203403 0.08305314  33.737847 1.611598e-249 nb      TRUE  26953.29908 16.478130
  # 4  chr15_61983466_61983605 chr15_61983466_61983605    383161     0.01515152  fixed      cond  (Intercept) 2.77292469 0.04600343  60.276471  0.000000e+00 nb      TRUE   2733.21045 16.005376
  # 5  chr15_61983326_61983465 chr15_61983466_61983605    383161     0.01515152  fixed      cond  (Intercept) 2.58336025 0.04323763  59.747964  0.000000e+00 nb      TRUE   5931.68101 13.241558
  # 6  chr15_61986956_61987105 chr15_61986956_61987105    383161     0.01515152  fixed      cond  (Intercept) 2.44032274 0.06643104  36.734675 2.042501e-295 nb      TRUE   2149.53659 11.476744
  #i=which(cell_ids=="chr15_61984436_61985515 chr15_61985676_61987145")
  #i=which(cell_ids=="chr15_61983326_61983465 chr15_61983326_61983465")
  i=which(cell_ids=="chr15_61986976_61987125 chr15_61987126_61988445")
  #browser()
    nb_mods<-app_fun(1:length(cell_ids),FUN=function(i) {
      #check if there are at least two obs per condition.
      if(lasso1d_filter){if(cell_ids[i] %in% dropped_covariates){return()}}
      tryCatch({
        #const used to avoid zero variance and failure to compare very small values below the winsorization threshold.
        #citation for this idea: Khan, J. A., Van Aelst, S., & Zamar, R. H. (2007). Robust Linear Model Selection Based on Least Angle Regression. Journal of the American Statistical Association, 102(480), 1289–1299. https://doi.org/10.1198/016214507000000950
        if(restricted_mnase){
        cell_df=long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt"),trtlr_cell_id==cell_ids[i])     %>% dplyr::mutate(value=ifelse(value<thresh,thresh,value),value=ifelse(value<=0,thresh,value))  %>% dplyr::mutate(value=as.integer(value)) %>% dplyr::mutate(value=ifelse(value>cap,cap,value))  
        } else{
        cell_df=long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt"))  %>% dplyr::mutate(value=ifelse(trtlr_cell_id==cell_ids[i],value,NA))   %>% dplyr::mutate(value=ifelse(value<thresh,thresh,value),value=ifelse(value<=0,thresh,value)) %>% dplyr::mutate(value=as.integer(value)) %>% dplyr::mutate(value=ifelse(value>cap,cap,value))
        }
        if(length(unique(cell_df$treatment))<2 ){
          message_parallel(paste0("less than two conditions for model",i," of ",length(cell_ids)))
          return();
        }
        if(min(table(droplevels(as.factor(cell_df$treatment))))<2){
          message_parallel(paste0("at least one condition has less than two observations for model",i," of ",length(cell_ids)))
          return()
        }
      #browser()
      gam_formula=base_formula
      
      if(offset){  gam_formula=paste0(gam_formula," + offset(log(total_reads))")}
      #gamTMB::gamTMB(value ~   treatment + offset(log(total_reads)) + splines::bs(mnase_bias),data= cell_df,family=fam,control=control_args)->nb_mod
         
      #else{
        #gamTMB::gamTMB(value ~   treatment ,data= cell_df,family=fam,control=control_args)->nb_mod
        #}
      if(mnase){
        gam_formula=paste0(gam_formula," + s(mnase,bs='ad')")
      }
      if(extra_formula_terms!=""){
        gam_formula=paste0(gam_formula," + ",extra_formula_terms)
      }
      #browser()
      # if(fam=="nb"){
      basepairC::message_parallel(paste0("model formula for model",i," of ",length(cell_ids),": ",gam_formula))
      if(ridge){
      mgcv::gam(as.formula(gam_formula),data= cell_df,family=eval(bquote(mgcv::nb(link = .(link)))),select=T,control=list(nthreads=single_model_cores),paraPen=list(treatment = list( diag(1))))->nb_mod } else{
        mgcv::gam(as.formula(gam_formula),data= cell_df,family=eval(bquote(mgcv::nb(link = .(link)))),select=T,control=list(nthreads=single_model_cores))->nb_mod #,weights = 1/log10(mnase)
      }
      
      
        # }
      # if(fam=="tw"){
      #   mgcv::gam(as.formula(gam_formula),data= cell_df,family=mgcv::tw(link=link),paraPen=list(treatment = list( diag(1))),select=T)->nb_mod
      # }
      if(bootstrap){
        nb_mod=bootstrap_gam_nb_hurdle(cell_name=cell_ids[i],formula_text=gam_formula,long_ctrl_trt_df,fam=gamTMB::truncated_nbinom2(),n_straps=1e3,n.cores=parallel::detectCores()/4,strap_size=NULL,return_full_results=F)
      }
      message_parallel(paste0("convergence; finished model",i," of ",length(cell_ids)))
      return(nb_mod)
      },error=function(e){warning(paste0(cell_ids[i]," failed to converge",i));return();})
      message_parallel(paste0("no convergence; finished model",i," of ",length(cell_ids)))
    })
    #browser()
    #if(return_models){return(list(nb_mods=nb_mods,zi_mods=zi_mods))}
    #browser()
    nb_mods_tidy_l<-app_fun(1:length(cell_ids),FUN=function(i) {
      if(is.null(nb_mods[[i]])){return(
        data.frame(term=c("(Intercept)","treatmenttrt"),estimate=c(0,0),std.error=c(1,1),statistic=c(0,0),p.value=c(1,1),model="nb",cell_id=cell_ids[i],AIC=0,converged=FALSE)
      )}
      if(any(class(nb_mods[[i]]) %in% "data.frame")){
        
        nb_mods[[i]] %>% dplyr::mutate(model="nb",cell_id=cell_ids[i],converged=T,AIC=0)->nb_mod_tidy
        
      }
      tryCatch({
        #browser()
      #broom.mixed::tidy(nb_mods[[i]])  %>% dplyr::mutate(model="nb",cell_id=cell_ids[i],converged=TRUE)-> nb_mod_tidy
      if(robust_cov){
      lmtest::coeftest(nb_mods[[i]],vcov. = mgcv::vcov.gam(nb_mods[[i]],sandwich=T,freq=T)) %>% broom::tidy() %>% dplyr::mutate(model="nb",cell_id=cell_ids[i])-> nb_mod_tidy} else{
        broom::tidy(nb_mods[[i]],parametric=T)  %>% dplyr::mutate(model="nb",cell_id=cell_ids[i])-> nb_mod_tidy
      }
        nb_mod_tidy$AIC=AIC(nb_mods[[i]])
        nb_mod_tidy$converged=T
      return(nb_mod_tidy)
      },error=function(e){warning(paste0(cell_ids[i]," failed to converge"));
        message_parallel(paste0(cell_ids[i]," failed to converge",i));
        return(
        data.frame(effect="fixed",component="cond",term=c("(Intercept)","treatmenttrt"),estimate=c(0,0),std.error=c(1,1),statistic=c(0,0),p.value=c(1,1),model="nb",cell_id=cell_ids[i],converged=FALSE,AIC=0)
        );})
    }) # %>% do.call(.,dplyr::bind_rows())  # %>% do.call(rbind,.) #%>% dplyr::group_by(cell_id) %>% dplyr::mutate(cell_area_frac=cell_area/sum(cell_area)) %>% dplyr::ungroup() #,mc.cores=min(length(cell_ids),parallel::detectCores()) #
    #%>% data.table::rbindlist()
    #browser()
    nb_mods_tidy=do.call(dplyr::bind_rows,nb_mods_tidy_l)
    dplyr::bind_rows(zi_mods_tidy,nb_mods_tidy) %>% dplyr::mutate(exp_est=exp(estimate),cell_area=get_cell_area(cell_id))   -> complete_model_results
    
    #lasso_results %>% dplyr::arrange(-exp_est)
    complete_model_results %>% dplyr::select(cell_id,cell_area) %>% dplyr::distinct() %>% dplyr::pull(cell_area) %>% sum()->total_cell_area
    complete_model_results %>% dplyr::mutate(cell_area_frac=cell_area/total_cell_area) %>% dplyr::select(cell_id,cell_area,cell_area_frac,dplyr::everything()) %>% dplyr::arrange(-exp_est) -> complete_model_results_with_cell_area
  #}
  # if(mc.cores>1){
  #   
  #   zi_mods<-app_fun(1:length(cell_ids),FUN=function(i) {
  #     basepaiC::message_parallel(paste0("zi mod",i," of ",length(cell_ids)))
  #     gamTMB::gamTMB(value==0 ~   treatment  ,data=long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt")  & trtlr_cell_id==cell_ids[i]) %>% dplyr::mutate(value=as.integer(value)) ,family=binomial)->zi_mod
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
  #     gamTMB::gamTMB(value ~   treatment  ,data= long_ctrl_trt_df_cell,family=fam)->nb_mod
  #     return(nb_mod)
  #     },error=function(e){warning(paste0(cell_ids[i]," failed to converge"));return();})
  #   },mc.cores=mc.cores)
  #   if(return_models){return(list(nb_mods=nb_mods,zi_mods=zi_mods))}
  #   nb_mods_tidy<-app_fun(1:length(cell_ids),FUN=function(i) {
  #     if(is.null(nb_mods[[i]])){return(
  #       data.frame(term=c("(Intercept)","treatmenttrt"),component="cond",estimate=c(0,0),std.error=c(1,1),statistic=c(0,0),p.value=c(1,1),converged=FALSE,model="nb",cell_id=cell_ids[i])
  #     )}
  #     tryCatch({
  #     broom.mixed::tidy(nb_mods[[i]])  %>% dplyr::mutate(model="nb",cell_id=cell_ids[i])-> nb_mod_tidy
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
  output_summary=complete_model_results_with_cell_area %>% tidyr::separate(cell_id,sep="_| ",into=c("chr2","start2","end2","chr1","start1","end1")) %>% dplyr::mutate(p.value=Rmpfr::pnorm(abs(statistic),lower.tail = F)) %>% dplyr::mutate(padj=p.adjust(p.value,method="fdr")) %>% dplyr::filter(model=="nb",term=="treatmenttrt") 

  plot(output_summary  %>% ggplot(aes(xmin=as.numeric(start1),xmax=as.numeric(end1),ymin=as.numeric(start2),ymax=as.numeric(end2),fill=ifelse(padj<0.05,(abs(statistic)),0))) + geom_rect() + scale_fill_viridis_c(direction = plot_direction,option="magma",oob=scales::squish)+ theme_minimal() + labs(title="truncated negative binomial gam Wald statistics (W=β/SE), by cell id",x="Start1",y="Start2")  ) + ggplot2::coord_equal() + ggplot2::theme(legend.position="bottom")
  
  plot(output_summary %>% dplyr::mutate(rate_ratio=as.numeric((exp(estimate)) ) ) %>% dplyr::arrange((rate_ratio))  %>% ggplot(aes(xmin=as.numeric(start1),xmax=as.numeric(end1),ymin=as.numeric(start2),ymax=as.numeric(end2),fill=rate_ratio)) + geom_rect() + scale_fill_viridis_c(option="magma",direction=plot_direction)+ theme_minimal() + labs(title="truncated negative binomial gam rate ratio values, by cell id",x="Start1",y="Start2")) + ggplot2::coord_equal() + ggplot2::theme(legend.position="bottom")
  #if(interactive){plotly::ggplotly(pval_ratio_plot) %>% plotly::toWebGL()} else {print(pval_ratio_plot)}
  output_summary=complete_model_results_with_cell_area %>% tidyr::separate(cell_id,sep="_| ",into=c("chr2","start2","end2","chr1","start1","end1"),remove = F)  %>% dplyr::mutate(p.value=Rmpfr::pnorm(abs(statistic),lower.tail = F)) %>% dplyr::mutate(padj=p.adjust(p.value,method="fdr")) %>% dplyr::arrange(term,-abs(statistic))
  #browser()
  if(return_models==F){
  return(list(output_summary=output_summary,cell_ids=cell_ids,complete_model_results_with_cell_area=complete_model_results_with_cell_area))} else{
    names(zi_mods)=cell_ids
    names(nb_mods)=cell_ids
    return(list(output_summary=output_summary,cell_ids=cell_ids,complete_model_results_with_cell_area=complete_model_results_with_cell_area,nb_mods=nb_mods,zi_mods=zi_mods))
  }
}
get_cell_area=function(cell_id) {
  (CNVScope::underscored_pos_to_GRanges(strsplit(cell_id," ")[[1]][[1]])@ranges@width)*(CNVScope::underscored_pos_to_GRanges(strsplit(cell_id," ")[[1]][[2]])@ranges@width) ->area
  return(area)
}