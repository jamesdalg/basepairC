#' Plot Submatrices for Different Treatments
#'
#' This function takes a dataframe of submatrices, filters by treatment types, and generates heatmaps
#' for each treatment type using `ggplot2` and `scattermore`. The plots are then combined and returned.
#' The function is specifically designed to visualize comparisons across control and treatment groups.
#'
#' @param submat_df A dataframe containing the submatrices with columns for start positions (`start1`, `start2`),
#'                  treatment types, and values to be plotted. The expected columns are `start1`, `start2`,
#'                  `treatment`, and `value`.
#'
#' @return A combined ggplot object containing the heatmaps for each treatment type.
#'
#' @importFrom ggplot2 ggplot aes labs scale_fill_viridis_c scale_color_viridis_c theme_minimal facet_wrap
#' @importFrom scattermore geom_scattermore
#'
#' @examples
#' # Assuming 'submat_data' is a dataframe with the appropriate structure:
#' submat_data <- data.frame(
#'   start1 = sample(1:100, 20),
#'   start2 = sample(1:100, 20),
#'   value = runif(20),
#'   treatment = rep(c("ctrl", "trt", "trtlr"), length.out = 20)
#' )
#' plot_submats(submat_data)
#'
plot_submats=function(submat_df) {
  submat_df %>% dplyr::filter(treatment=="ctrl")  %>% ggplot(aes(x=start1,y=start2,fill=value,color=value)) + scattermore::geom_scattermore() + scale_fill_viridis_c(option = "inferno") + scale_color_viridis_c(option = "inferno") + theme_minimal() + labs(title="Original data heatmap",x="Start1",y="Start2") + facet_wrap(~treatment) ->ctrl_hm
  submat_df %>% dplyr::filter(treatment=="trt")  %>% ggplot(aes(x=start1,y=start2,fill=value,color=value)) + scattermore::geom_scattermore() + scale_fill_viridis_c(option = "inferno") + scale_color_viridis_c(option = "inferno") + theme_minimal() + labs(title="Original data heatmap",x="Start1",y="Start2") + facet_wrap(~treatment) ->trt_hm
  submat_df %>% dplyr::filter(treatment=="trtlr")  %>% ggplot(aes(x=start1,y=start2,fill=value,color=value)) + scattermore::geom_scattermore() + scale_fill_viridis_c(option = "inferno") + scale_color_viridis_c(option = "inferno") + theme_minimal() + labs(title="Original data heatmap",x="Start1",y="Start2") + facet_wrap(~treatment)->trtlr_hm
  return(ctrl_hm + trt_hm + trtlr_hm)
}
