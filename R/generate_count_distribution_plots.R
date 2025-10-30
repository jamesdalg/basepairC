#' Generate Discrete Distribution Plots for Control and Treatment Groups
#'
#' This function applies the `descdist` function from the `fitdistrplus` package to generate
#' discrete distribution plots for two subsets of data: control and treatment groups.
#' Each group's data is filtered from the main dataframe, and the distribution of the
#' `value` column is analyzed.
#'
#' @param df A dataframe containing at least the columns `treatment` and `value`.
#'           The `treatment` column should have at least two levels: "ctrl" for control
#'           and "trt" for treatment.
#'
#' @return Returns a list containing the distribution plot objects for both the control
#'         (control_plot) and treatment (treatment_plot) groups. Each object is the result
#'         of the `descdist` function, providing a comprehensive analysis of the empirical
#'         distribution for each group.
#'
#' @examples
#' # Assuming `data` is a dataframe with `treatment` and `value` columns
#' plots <- generate_count_distribution_plots_discrete(data)
#' # You can then access each plot object via plots$control_plot and plots$treatment_plot
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @import fitdistrplus
#' @export
generate_count_distribution_plots_discrete=function(df){
  control_plot <- fitdistrplus::descdist(df %>% dplyr::filter(treatment=="ctrl") %>% dplyr::select(value) %>% unlist(),boot=100,discrete=T)
  treatment_plot <- fitdistrplus::descdist(df %>% dplyr::filter(treatment=="trt") %>% dplyr::select(value) %>% unlist(),boot=100,discrete=T)
  return(list(control_plot=control_plot, treatment_plot=treatment_plot))
}
