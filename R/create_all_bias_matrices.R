#' Create Bias Matrices for Multiple Columns
#'
#' This function generates a list of bias matrices for each column in the input data table, excluding the specified columns: `"chr"`, `"start"`, `"end"`, and `"seq"`. The bias matrices are created by calling the `create_bias_matrix_from_column` function for each column.
#'
#' @param bias_dt A data.table containing the bias data. It must include columns named `"chr"`, `"start"`, and `"end"`, and optionally a column `"seq"`. All other columns will be used to generate the bias matrices.
#' @return A named list of bias matrices, where each matrix corresponds to one of the columns in `bias_dt` that is not `"chr"`, `"start"`, `"end"`, or `"seq"`. The names of the list will correspond to the column names in `bias_dt`.
#' @details 
#' The function iterates over the columns of `bias_dt` excluding `"chr"`, `"start"`, `"end"`, and `"seq"`. For each remaining column, it applies the `create_bias_matrix_from_column` function, which constructs the bias matrix for that column. 
#' 
#' @seealso \code{\link{create_bias_matrix_from_column}} for details on how the bias matrices are generated.
#' @importFrom pbapply pblapply
#' 
#' @examples
#' \dontrun{
#'   bias_dt <- data.table::data.table(chr = rep("chr1", 100),
#'                                     start = seq(1, 1000, by=10),
#'                                     end = seq(10, 1000, by=10),
#'                                     bias1 = rnorm(100),
#'                                     bias2 = runif(100))
#'   bias_matrices <- create_all_bias_matrices(bias_dt)
#'   print(names(bias_matrices))
#' }
#' @export
create_all_bias_matrices=function(bias_dt){
  bias_cols=setdiff(names(bias_dt),c("chr","start","end","seq"))
  bias_matrices=pbapply::pblapply(bias_cols,function(col_name){
    create_bias_matrix_from_column(bias_dt[[col_name]],bias_dt$chr,bias_dt$start,bias_dt$end,minmaxnorm=F)
  })
  names(bias_matrices)=bias_cols
  return(bias_matrices)
}
