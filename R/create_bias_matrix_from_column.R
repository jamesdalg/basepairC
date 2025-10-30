#' Create a Bias Matrix from a Column of Bias Values
#'
#' This function generates a bias matrix from a vector of bias values, with optional min-max normalization. The matrix is created by calculating the outer product of the bias values, and the resulting matrix has row and column names constructed from the chromosome, start, and end values.
#'
#' @param bias_values A numeric vector of bias values used to construct the bias matrix.
#' @param chr_vals A vector of chromosome identifiers corresponding to the bias values.
#' @param start_vals A vector of start positions corresponding to the bias values.
#' @param end_vals A vector of end positions corresponding to the bias values.
#' @param minmaxnorm Logical. If `TRUE`, the bias values are normalized using min-max normalization before constructing the bias matrix. Default is `FALSE`.
#' @return A square bias matrix where the dimensions correspond to the length of `bias_values`. The row and column names are derived from the combination of `chr_vals`, `start_vals`, and `end_vals`.
#' @details 
#' The function first checks whether `minmaxnorm` is `TRUE`. If so, the bias values are normalized using the `basepairC::minmaxnorm` function. Then, a bias row is constructed by replicating the (optionally normalized) bias values. This row is transposed to form a bias column, and the outer product of the row and column is used to create the bias matrix. The row and column names of the matrix are set using the `chr_vals`, `start_vals`, and `end_vals`.
#' 
#' @importFrom magrittr %>%
#' @examples
#' 
#' \dontrun{
#'   bias_values <- rnorm(100)
#'   chr_vals <- rep("chr1", 100)
#'   start_vals <- seq(1, 1000, by=10)
#'   end_vals <- seq(10, 1000, by=10)
#'   bias_matrix <- create_bias_matrix_from_column(bias_values, chr_vals, start_vals, end_vals, minmaxnorm=TRUE)
#'   print(dim(bias_matrix))
#' }

#' @export
create_bias_matrix_from_column=function(bias_values,chr_vals,start_vals,end_vals,minmaxnorm=F){
  if(minmaxnorm){
    matrix(rep(bias_values %>% basepairC::minmaxnorm(),length(bias_values)),nrow=length(bias_values),byrow=T) ->bias_row} else{
      matrix(rep(bias_values,length(bias_values)),nrow=length(bias_values),byrow=T) ->bias_row
    }
  bias_row %>% t()->bias_col
  bias_matrix=(bias_row*bias_col)
  rownames(bias_matrix)<-colnames(bias_matrix)<-paste0(chr_vals,"_",start_vals,"_",end_vals)
  return(bias_matrix)
}
