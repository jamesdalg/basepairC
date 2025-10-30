#' @title gam2d_negbin_surface_model
#' @description This function fits a 2D GAM model to the control/trt ratios using a negative binomial distribution.
#' @param long_ctrl_trt_df A long-format data frame with control/treatment ratio data.
#' @param knots The number of knots to use in the 2D GAM model.
#' @param plot A logical indicating whether to plot the 2D GAM model.
#' @importFrom mgcv gam te predict.gam
#' @importFrom ggplot2 geom_raster scale_color_viridis_c theme_minimal labs stat_contour scale_fill_viridis_c geom_contour aes ggplot
#' @importFrom magrittr `%>%`
#' @return pred_gam2d A data frame with the predicted values and standard errors from the 2D GAM model.
#' @export gam2d_negbin_surface_model
gam2d_negbin_surface_model=function(long_ctrl_trt_df,knots=100,plot=T){
  cell_ids<-long_ctrl_trt_df %>% dplyr::select(trtlr_cell_id) %>% unique() %>% dplyr::pull()
  gam_fun=value~trtlr_cell_id+te(start1,start2,bs=c("ts"),k=knots,d=2)
  gam2d<-gam(gam_fun,data=long_ctrl_trt_df %>% dplyr::filter(treatment=="trtlr"),family=nb(link="log"),nthreads=parallel::detectCores())
  predict.gam(gam2d,newdata=long_ctrl_trt_df %>% dplyr::filter(treatment=="trtlr"),type="response",se.fit = T) -> pred_gam2d
  if(plot){
    print(long_ctrl_trt_df %>% dplyr::filter(treatment=="trtlr") %>% dplyr::mutate(se_gam2d=pred_gam2d$se.fit,pred_gam2d=pred_gam2d$fit) %>% ggplot(aes(x=start1,y=start2,fill=pred_gam2d/se_gam2d,z=pred_gam2d/se_gam2d)) + geom_raster() + scale_color_viridis_c(option = "inferno") + theme_minimal() + labs(title="GAM plot for trtlr cell id",x="Start1",y="Start2") + stat_contour()  + scale_fill_viridis_c(option = "inferno") + geom_contour(aes(z=pred_gam2d/se_gam2d)))->gam_plot
  }
  return(pred_gam2d)
}