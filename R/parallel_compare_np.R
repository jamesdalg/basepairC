#' @title parallel_compare (non-parametric)
#' @description This function performs a parallel comparison of the treatment effect of a control and treatment group with Wilcoxon & permutation tests.
#' @param long_ctrl_trt_df A long-format data frame with control and treatment ratio data.
#' @param mc.cores The number of CPU cores to use for parallel processing (default is 1).
#' @param interactive A logical value indicating whether to plot the results interactively (default is FALSE).
#' @return A data frame with the treatment effect comparison results.
#' @importFrom magrittr `%>%`
#' @import patchwork
#' @examples
#' dontrun{
#' parallel_compare(long_ctrl_trt_df, mc.cores = 2, interactive = TRUE) 
#' }
#' @export
parallel_compare_np=function(long_ctrl_trt_df,mc.cores=1,interactive=F,alternative="less"){
  cell_ids=long_ctrl_trt_df %>% dplyr::filter(treatment=="trtlr") %>% dplyr::select(trtlr_cell_id) %>% unique() %>% dplyr::pull()
  pbapply_funs=list(pblapply=pbapply::pblapply, pbmcapply=pbmcapply::pbmclapply)
  if(mc.cores==1){
    app_fun=pbapply_funs$pblapply
  }    else {
    app_fun=pbapply_funs$pbmcapply
  }
    np_tidy<-app_fun(1:length(cell_ids),FUN=function(i) {

      message_parallel("Processing cell id ",cell_ids[i])
      message_parallel("Processing cell ",i," of ",length(cell_ids))
      cell_df=long_ctrl_trt_df %>% dplyr::filter(treatment %in% c("ctrl","trt"),trtlr_cell_id==cell_ids[i]) 
      mean_diff =c((cell_df %>% dplyr::filter(treatment=="trt") %>% dplyr::pull(value))-(cell_df %>% dplyr::filter(treatment=="ctrl") %>% dplyr::pull(value))) %>% mean()
      permTS_tidy=data.frame(estimate=mean_diff,statistic=NA,p.value=NA,model="permTS",cell_id=cell_ids[i],alternative=paste0(alternative))
      wilcox_tidy=data.frame(statistic=NA,p.value=NA,estimate=mean_diff,model="wilcoxon",cell_id=cell_ids[i],alternative=paste0(alternative))
      tryCatch({
      #statistic p.value method                                               alternative model  cell_id                                         estimate
      wilcox.test(cell_df %>% dplyr::filter(treatment=="trt") %>% dplyr::pull(value),cell_df %>% dplyr::filter(treatment=="ctrl") %>% dplyr::pull(value),alternative=alternative,paired=T)->wilcox_out
      wilcox_out %>% broom::tidy() %>% dplyr::mutate(model="wilcoxon",cell_id=cell_ids[i],alternative=paste0(alternative),estimate=mean_diff) -> wilcox_tidy
      },error=function(e){message_parallel("Error in cell id ",cell_ids[i],": ",e$message);
        wilcox_tidy=data.frame(statistic=NA,p.value=NA,estimate=mean_diff,model="wilcoxon",cell_id=cell_ids[i],alternative=paste0(alternative))})
      
      tryCatch({
      #estimate statistic p.value method                                          alternative model  cell_id 
      perm::permTS(x=cell_df %>% dplyr::filter(treatment=="trt") %>% dplyr::pull(value),y=cell_df %>% dplyr::filter(treatment=="ctrl") %>% dplyr::pull(value),alternative=alternative)->permTS_out
      permTS_out %>% broom::tidy() %>% dplyr::mutate(model="permTS",cell_id=cell_ids[i],alternative=paste0(alternative)) -> permTS_tidy
      },error=function(e){message_parallel("Error in cell id ",cell_ids[i],": ",e$message);
        permTS_tidy=data.frame(estimate=mean_diff,statistic=NA,p.value=NA,model="permTS",cell_id=cell_ids[i],alternative=paste0(alternative))})
      #do the same for ks.test()
      
      ks_tidy=data.frame(statistic=NA,p.value=NA,estimate=mean_diff,model="ks",cell_id=cell_ids[i],alternative=paste0(alternative))
      tryCatch({
      ks.test(cell_df %>% dplyr::filter(treatment=="trt") %>% dplyr::pull(value),cell_df %>% dplyr::filter(treatment=="ctrl") %>% dplyr::pull(value),alternative=alternative)->ks_out
      ks_out %>% broom::tidy() %>% dplyr::mutate(estimate=mean_diff,model="ks",cell_id=cell_ids[i]) ->ks_tidy
      },error=function(e){message_parallel("Error in cell id ",cell_ids[i],": ",e$message);
        ks_tidy=data.frame(statistic=NA,p.value=NA,estimate=mean_diff,model="ks",cell_id=cell_ids[i],alternative=paste0(alternative))})
      
      
      
      return(dplyr::bind_rows(permTS_tidy,wilcox_tidy,ks_tidy))
      #return(permTS_tidy)
      #combine permutation and wilcoxon outputs into tidy form.
      #glmmTMB::glmmTMB(value ~   condition  ,data=long_ctrl_trt_df %>% dplyr::filter(condition %in% c("ctrl","trt") & treatment=="lasso" & trtlr_cell_id==cell_ids[i] & value>thresh) %>% dplyr::mutate(value=as.integer(value)) ,family=glmmTMB::truncated_nbinom2())->nb_mod

    }) %>% do.call(dplyr::bind_rows,.) #data.table::rbindlist() #this works better because dplyr is slower but matches the columns when there is a different order possible (e.g. when one test converges and the other does not).
    np_tidy %>% dplyr::mutate(exp_est=exp(estimate),cell_area=get_cell_area(cell_id))   -> complete_model_results
    
    #lasso_results %>% dplyr::arrange(-exp_est)
    complete_model_results %>% dplyr::select(cell_id,cell_area) %>% dplyr::distinct() %>% dplyr::pull(cell_area) %>% sum()->total_cell_area
    complete_model_results %>% dplyr::mutate(cell_area_frac=cell_area/total_cell_area) %>% dplyr::select(cell_id,cell_area,cell_area_frac,dplyr::everything()) %>% dplyr::arrange(-exp_est) -> complete_model_results_with_cell_area
    # complete_model_results_with_cell_area
    #make a mean diff plot
    #works
    
    
    #expermental (works but only useful for interactivity really)
    # complete_model_results_with_cell_area %>% tidyr::pivot_wider(names_from=model,values_from =p.value) %>% tidyr::unnest() %>% replace_na(replace=list(permTS=1,wilcoxon=1))    %>% dplyr::select(cell_id,estimate,wilcoxon,permTS) %>% dplyr::distinct() %>% tidyr::separate(col = cell_id,sep = " |_",into=c('chr1','start1','end1','chr2',"start2","end2"),remove = F) %>%
    #   dplyr::mutate(start1=as.numeric(start1),start2=as.numeric(start2),end1=as.numeric(end1),end2=as.numeric(end2)) %>%
    #   dplyr::filter(start1<=start2)->plot_df
    # long_ctrl_trt_df %>% dplyr::inner_join(plot_df,by=c("trtlr_cell_id"="cell_id")) %>% dplyr::distinct() %>% dplyr::rename(cell_id=trtlr_cell_id,start1=start1.x,start2=start2.x)  %>%
    #   #tidyr::separate(col = cell_id,sep = " |_",into=c('chr1','start1','end1','chr2',"start2","end2"),remove = F) %>%
    #   dplyr::select(cell_id,estimate,start1,start2,end2,wilcoxon,permTS) %>% dplyr::distinct()  %>%
    #   dplyr::mutate(start1=as.numeric(start1),start2=as.numeric(start2)) %>%
    #   dplyr::filter(start1<=start2) -> plot_df_merged
    # plot_df_merged %>% ggplot2::ggplot(ggplot2::aes(x=start1,y=start2,color=estimate,fill=estimate,text=paste0("wilcoxon p-value=",wilcoxon,"/n permutation test p-value=",permTS))) +
    #   geom_raster() + scale_color_viridis_c() + scale_fill_viridis_c() + theme_minimal() + 
    #   labs(title="Mean diff results for trtlr cell id",x="Start1",y="Start2") -> mean_diff_plot_merged
    
    
    #works (below)
    complete_model_results_with_cell_area %>% dplyr::select(cell_id,estimate) %>% dplyr::distinct() %>% tidyr::separate(col = cell_id,sep = " |_",into=c('chr1','start1','end1','chr2',"start2","end2"),remove = F) %>% 
      dplyr::mutate(start1=as.numeric(start1),start2=as.numeric(start2),end1=as.numeric(end1),end2=as.numeric(end2)) %>% 
      dplyr::filter(start1<=start2)->plot_df
    
    long_ctrl_trt_df %>% dplyr::inner_join(plot_df %>% dplyr::rename(cell_start1=start1,cell_start2=start2,cell_end1=end1,cell_end2=end2),by=c("trtlr_cell_id"="cell_id")) -> plot_df_merged
    plot_df %>% ggplot2::ggplot(ggplot2::aes(xmin=start1,xmax=end1,ymin=start2,ymax=end2,color=estimate,fill=estimate)) +
      ggplot2::geom_rect() + ggplot2::scale_color_viridis_c() + ggplot2::scale_fill_viridis_c() + ggplot2::theme_minimal() + 
      ggplot2::labs(title="Mean diff results for trtlr cell id",x="Start1",y="Start2") -> mean_diff_plot
    plot_df_merged %>% ggplot2::ggplot(ggplot2::aes(x=start1,y=start2,color=estimate,fill=estimate)) +  #,text=paste0("wilcoxon p-value=",wilcoxon,"/n permutation test p-value=",permTS)
      ggplot2::geom_raster() + ggplot2::scale_color_viridis_c() + ggplot2::scale_fill_viridis_c() + ggplot2::theme_minimal() + 
      ggplot2::labs(title="Mean diff results for trtlr cell id",x="Start1",y="Start2") -> mean_diff_plot_merged
    print(mean_diff_plot_merged)
    #complete_model_results_with_cell_area %>% dplyr::select(-method) %>% tidyr::pivot_wider(names_from=c("model"),values_from=c("statistic","p.value","alternative","exp_est")) %>% tidyr::unnest() 
    #complete_model_results_with_cell_area  %>% tidyr::separate(col = cell_id,sep = " |_",into=c('chr1','start1','end1','chr2',"start2","end2"),remove = F) %>% dplyr::filter(start1<start2,model=="wilcoxon") %>% dplyr::mutate(start1=as.numeric(start1),end1=as.numeric(end1),start2=as.numeric(end2))  %>% 
     # ggplot2::ggplot(aes(xmin=start1,xmax=start2,ymin=start1,ymax=start2,color=estimate,fill=estimate)) + geom_rect() + scale_color_viridis_c() + scale_fill_viridis_c() + theme_minimal() + labs(title="Wilcoxon results for trtlr cell id",x="Start1",y="Start2") -> wilcox_plot
    #complete_model_results_with_cell_area %>% dplyr::select(-method) %>% tidyr::pivot_wider(names_from=c("model"),values_from=c("statistic","p.value","alternative","exp_est")) %>% tidyr::unnest() 
   
    
  #   complete_model_results_with_cell_area %>% dplyr::filter(start1<=start2)  %>% dplyr::mutate(plotted_value=ifelse(p.value!=0,-log(p.value),-log(min(complete_merged_glm_results$p.value[complete_merged_glm_results$p.value>0]))) %>% CNVScope::signedRescale())
  # ) %>% ggplot(aes(x=start1,y=start2,color=plotted_value,fill=plotted_value)) + geom_raster() + scale_color_viridis_c() + scale_fill_viridis_c() + theme_minimal() + labs(title="GLM results for trtlr cell id",x="Start1",y="Start2") -> pval_ratio_plot
  # if(interactive){plotly::ggplotly(pval_ratio_plot) %>% plotly::toWebGL()} else {print(pval_ratio_plot)}
  output_df=complete_model_results_with_cell_area %>% dplyr::mutate(padj_overall=p.adjust(p.value,method="fdr")) %>% dplyr::group_by(model) %>% dplyr::mutate(padj_model=p.adjust(p.value,method="fdr")) %>% dplyr::arrange(padj_overall,-estimate) %>% tidyr::unnest(cols=c())
  #p-value plot
  complete_model_results_with_cell_area  %>% dplyr::filter(model=="ks") %>% dplyr::select(cell_id,estimate,p.value) %>% dplyr::distinct() %>% tidyr::separate(col = cell_id,sep = " |_",into=c('chr1','start1','end1','chr2',"start2","end2"),remove = F) %>% 
    dplyr::mutate(start1=as.numeric(start1),start2=as.numeric(start2),end1=as.numeric(end1),end2=as.numeric(end2)) %>% 
    dplyr::filter(start1<=start2)->plot_df_ks
  long_ctrl_trt_df %>% dplyr::inner_join(plot_df_ks %>% dplyr::rename(cell_start1=start1,cell_start2=start2,cell_end1=end1,cell_end2=end2),by=c("trtlr_cell_id"="cell_id")) %>% dplyr::mutate(padj=p.adjust(p.value),sig=as.integer(padj<0.05))-> plot_df_merged_ks
  plot_df_merged_ks  %>% ggplot2::ggplot(ggplot2::aes(x=start1,y=start2,color=sig,fill=sig)) +  #,text=paste0("wilcoxon p-value=",wilcoxon,"/n permutation test p-value=",permTS)
    ggplot2::geom_raster() + ggplot2::scale_color_viridis_c() + ggplot2::scale_fill_viridis_c() + ggplot2::theme_minimal() + 
    ggplot2::labs(title="significant 2D segments using Kolmogorov-Smirnov test (nonparametric)",x="Start1",y="Start2") + ggplot2::theme(legend.position="bottom") + ggplot2::coord_equal()-> padj_plot_ks
  #experimental
  

  #plot(padj_plot_ks)
  #do the same for wilcoxon
  complete_model_results_with_cell_area  %>% dplyr::filter(model=="wilcoxon") %>% dplyr::select(cell_id,estimate,p.value) %>% dplyr::distinct() %>% tidyr::separate(col = cell_id,sep = " |_",into=c('chr1','start1','end1','chr2',"start2","end2"),remove = F) %>% 
    dplyr::mutate(start1=as.numeric(start1),start2=as.numeric(start2),end1=as.numeric(end1),end2=as.numeric(end2)) %>% 
    dplyr::filter(start1<=start2)->plot_df_wilcox
  long_ctrl_trt_df %>% dplyr::inner_join(plot_df_wilcox %>% dplyr::rename(cell_start1=start1,cell_start2=start2,cell_end1=end1,cell_end2=end2),by=c("trtlr_cell_id"="cell_id")) %>% dplyr::mutate(padj=p.adjust(p.value),sig=as.integer(padj<0.05))-> plot_df_merged_wilcox
  plot_df_merged_wilcox  %>% ggplot2::ggplot(ggplot2::aes(x=start1,y=start2,color=sig,fill=sig)) +  #,text=paste0("wilcoxon p-value=",wilcoxon,"/n permutation test p-value=",permTS)
    ggplot2::geom_raster() + ggplot2::scale_color_viridis_c() + ggplot2::scale_fill_viridis_c() + ggplot2::theme_minimal() + 
    ggplot2::labs(title="significant 2D segments using paired wilcoxon test (nonparametric)",x="Start1",y="Start2") + ggplot2::theme(legend.position="bottom") + ggplot2::coord_equal()-> padj_plot_wilcox
  #do the same for permTS
  complete_model_results_with_cell_area  %>% dplyr::filter(model=="permTS") %>% dplyr::select(cell_id,estimate,p.value) %>% dplyr::distinct() %>% tidyr::separate(col = cell_id,sep = " |_",into=c('chr1','start1','end1','chr2',"start2","end2"),remove = F) %>% 
    dplyr::mutate(start1=as.numeric(start1),start2=as.numeric(start2),end1=as.numeric(end1),end2=as.numeric(end2)) %>% 
    dplyr::filter(start1<=start2)->plot_df_permTS
  long_ctrl_trt_df %>% dplyr::inner_join(plot_df_permTS %>% dplyr::rename(cell_start1=start1,cell_start2=start2,cell_end1=end1,cell_end2=end2),by=c("trtlr_cell_id"="cell_id")) %>% dplyr::mutate(padj=p.adjust(p.value),sig=as.integer(padj<0.05))-> plot_df_merged_permTS
  plot_df_merged_permTS  %>% ggplot2::ggplot(ggplot2::aes(x=start1,y=start2,color=sig,fill=sig)) +  #,text=paste0("wilcoxon p-value=",wilcoxon,"/n permutation test p-value=",permTS)
    ggplot2::geom_raster() + ggplot2::scale_color_viridis_c() + ggplot2::scale_fill_viridis_c() + ggplot2::theme_minimal() + 
    ggplot2::labs(title="significant 2D segments using permutation test (nonparametric)",x="Start1",y="Start2") + ggplot2::theme(legend.position="bottom") + ggplot2::coord_equal()-> padj_plot_permTS
  #plot(padj_plot_wilcox)
  #plot(padj_plot_permTS)
  #plot(padj_plot_ks)
  mean_diff_plot_merged / padj_plot_wilcox + padj_plot_permTS + padj_plot_ks -> output_plot
  plot(output_plot)
  return(list(output_df=output_df,mean_diff_plot=mean_diff_plot_merged,padj_plot_wilcox=padj_plot_wilcox,padj_plot_permTS=padj_plot_permTS,padj_plot_ks=padj_plot_ks,output_plot=output_plot))
}
#tidyr::separate(col = cell_id,sep = " |_",into=c('chr1','start1','end1','chr2',"start2","end2"),remove = F) %>%
# long_ctrl_trt_df %>% dplyr::inner_join(plot_df,by=c("trtlr_cell_id"="cell_id")) %>% 
#   dplyr::rename(start1=start1.x,start2=start2.x) %>% dplyr::select(trtlr_cell_id,estimate,start1,end1,start2,end2) 
# %>% dplyr::distinct() %>% dplyr::rename(cell_id=trtlr_cell_id) %>%
#   dplyr::mutate(start1=as.numeric(start1),start2=as.numeric(start2),end1=as.numeric(end1),end2=as.numeric(end2)) %>%
#   dplyr::filter(start1<=start2) -> plot_df_merged