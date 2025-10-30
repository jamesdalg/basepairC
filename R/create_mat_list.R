#' Create a List of Matrices from rr_output
#'
#' This function processes the output of `rr_output` and returns a list of matrices,
#' one for each sample in the input matrix. Each matrix will have specified row and 
#' column names and attributes, including `condition` and `sample_name`.
#'
#' @param rr_output A list containing at least a `sample_mat` matrix, with samples 
#'   in columns and observations in rows. Additionally, `underscored_positions_col` 
#'   and `underscored_positions_row` should contain the column and row names, respectively.
#' @param ncol_output Optional. The number of columns for the output matrices. If 
#'   not provided, the number of columns will be set to the square root of the 
#'   number of rows in `rr_output$sample_mat`.
#' @param nrow_output Optional. The number of rows for the output matrices. If 
#'   not provided, the number of rows will be set to the square root of the 
#'   number of rows in `rr_output$sample_mat`.
#' @param ncores Optional. The number of cores to use for parallel processing. 
#'   If set to "auto", the number of cores will be automatically detected based 
#'   on the available cores on the machine. Default is 1.
#'
#' @return A list of matrices. Each matrix corresponds to one sample from 
#'   `rr_output$sample_mat` and has dimensions `nrow_output` by `ncol_output`. 
#'   The matrices will also have the attributes `condition` and `sample_name`.
#'
#' @details Creates a list of matrices, with condition and sample name attribute information.
#'    The column and row names of each matrix are taken from 
#'   `rr_output$underscored_positions_col` and `rr_output$underscored_positions_row`, 
#'   respectively. Each matrix is stored with attributes for `condition` and 
#'   `sample_name`.
#'
#' @importFrom pbmcapply pbmclapply
#' @importFrom parallel detectCores
#' @export
#'
#' @examples
#' \dontrun{
#' rr_output <- list(
#'   sample_mat = matrix(rnorm(100), 10, 10), 
#'   underscored_positions_col = paste0("col_", 1:10), 
#'   underscored_positions_row = paste0("row_", 1:10)
#' )
#' create_mat_list(rr_output, ncores = "auto")
#' }

create_mat_list=function(rr_output,ncol_output=NULL,nrow_output=NULL,ncores=1){
  if(is.null(ncol_output) | is.null(nrow_output)){
    ncol_output<-nrow_output<-sqrt(nrow(rr_output$sample_mat))}
  if(ncores=="auto"){options(mc.cores = min(parallel::detectCores(), ncol(rr_output$sample_mat)))} else {options(mc.cores = ncores)}
  pbmcapply::pbmclapply(1:ncol(rr_output$sample_mat), function(i) {
    matrix(as.data.frame(rr_output$sample_mat)[,i], nr = nrow_output, nc = ncol_output) -> output_mat
    colnames(output_mat) <- rr_output$underscored_positions_col
    rownames(output_mat) <- rr_output$underscored_positions_row
    attributes(output_mat)$condition=rr_output$conditions[i]
    attributes(output_mat)$sample_name=colnames(rr_output$sample_mat)[i]
    return(output_mat)
  }) -> sample_mat_list
  return(sample_mat_list)
}