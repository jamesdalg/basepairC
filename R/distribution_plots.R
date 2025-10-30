#' @title Generate distribution plots
#' @description This function generates distribution plots for ctrl/trt ratios and their transformations.
#' @param df A long-format data frame with control and treatment ratio data.
#' @param treatment A character string indicating the treatment group.
#' @return output_list A list of distribution plots and two data frames with distribution statistics for before the transformation and after the transformation.
#' @importFrom magrittr `%>%`
#' @export
generate_distibution_plots=function(df,treatment){
  #distribution plots
  (df %>% dplyr::filter(treatment=="trtlr") %>% ggplot2::ggplot(aes(x=sqrt(value))) + ggplot2::geom_histogram(bins=100) + ggplot2::labs(title="Histogram of sqrt trtlr values",x="Value",y="Frequency"))
  (df %>% dplyr::filter(treatment=="trtlr") %>% ggplot2::ggplot(aes(x=sqrt(value))) + ggplot2::geom_histogram(bins=100) + ggplot2::labs(title="Histogram of trtlr values",x="Value",y="Frequency"))
  fitdistrplus::descdist(df %>% dplyr::filter(treatment=="trtlr") %>% dplyr::select(value) %>% unlist(),boot=11) ->distplot
  fitdistrplus::descdist(df %>% dplyr::filter(treatment=="trtlr") %>% dplyr::select(value) %>% unlist(),boot=100,discrete=T) ->distplot_discrete
  fitdistrplus::descdist(sqrt(df %>% dplyr::filter(treatment==treatment) %>% dplyr::pull(value) ),boot=100) ->distplot_log
  fitdistrplus::descdist(log(df %>% dplyr::filter(treatment==treatment) %>% dplyr::pull(value) +1),boot=100) ->distplot_log_plus_one
  fitdistrplus::descdist((df %>% dplyr::filter(treatment==treatment) %>% dplyr::pull(value))^2,boot=100) ->distplot_square
  fitdistrplus::descdist((df %>% dplyr::filter(treatment==treatment) %>% dplyr::pull(value))^(1/2),boot=100) ->distplot_sqrt
  df %>% dplyr::filter(treatment=="trtlr") %>% dplyr::group_by(trtlr_cell_id) %>% dplyr::summarize(mean=mean(value),var=var(value)) ->cell_means_vars
  df %>% dplyr::filter(treatment=="trtlr") %>% dplyr::group_by(trtlr_cell_id) %>% dplyr::summarize(mean=mean(2^value-1),var=var(2^value-1)) ->cell_means_before_transform
  cell_means_vars %>%  knitr::kable()
  cell_means_before_transform %>% knitr::kable()
  output_list=list(distplot,distplot_discrete,distplot_log,distplot_log_plus_one,distplot_square,distplot_sqrt,cell_means_vars,cell_means_before_transform)
  return(output_list)
}