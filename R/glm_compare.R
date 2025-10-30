#' @title parallel_compare
#' @description This function performs a parallel comparison of the treatment effect of a control and treatment group.
#' @param long_ctrl_trt_df A long-format data frame with control and treatment ratio data.
#' @param interactive A logical value indicating whether to plot the results interactively (default is FALSE).
#' @return A data frame with the treatment effect comparison results.
#' @importFrom magrittr `%>%`
#' @examples
#' dontrun{
#' glm_compare(long_ctrl_trt_df, interactive = TRUE)
#' }
#' @export
glm_compare=function(long_ctrl_trt_df,interactive=F) {
  glm(data=long_ctrl_trt_df %>% dplyr::filter(treatment=="trtlr"),formula=value ~ 0 + trtlr_cell_id, family=gaussian(link="identity")) ->glm_trtlr
  glm_trtlr %>% broom::tidy() %>% dplyr::filter(stringr::str_detect(term,"trtlr_cell_id")) %>% dplyr::mutate(trtlr_cell_id=stringr::str_replace(term,pattern="trtlr_cell_id",replacement="")) %>% dplyr::mutate(trt_wt_ratio=2^estimate-1,wt_trt_ratio=1/trt_wt_ratio) %>% tidyr::drop_na() %>% dplyr::inner_join(long_ctrl_trt_df) %>%  dplyr::mutate(plotted_val=ifelse(start1>start2,2^(estimate)-1,ifelse(p.value!=0,-log10(p.value),-log(0+.Machine$double.eps))))->merged_glm_results
  
  (
    merged_glm_results %>% dplyr::filter(start1>start2)  %>% dplyr::mutate(plotted_value=CNVScope::signedRescale((trt_wt_ratio^(1/4)))) 
  ) %>% dplyr::bind_rows(
    
    merged_glm_results %>% dplyr::filter(start1<=start2)  %>% dplyr::mutate(plotted_value=ifelse(p.value!=0,-log(p.value),-log(min(merged_glm_results$p.value[merged_glm_results$p.value>0]))) %>% CNVScope::signedRescale())
  ) %>% ggplot(aes(x=start1,y=start2,color=plotted_value,fill=plotted_value)) + geom_raster() + scale_color_viridis_c() + scale_fill_viridis_c() + theme_minimal() + labs(title="GLM results for trtlr cell id",x="Start1",y="Start2") -> pval_ratio_plot
  if(interactive){plotly::ggplotly(pval_ratio_plot) %>% plotly::toWebGL()} else {print(pval_ratio_plot)}
  return(merged_glm_results %>% dplyr::select(trtlr_cell_id,trt_wt_ratio,estimate,std.error,statistic,p.value) %>% dplyr::distinct() %>% dplyr::mutate(padj=p.adjust(Rmpfr::pnorm(estimate),method="fdr")))
}
