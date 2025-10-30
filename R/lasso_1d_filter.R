#' LASSO-Based Filtering of Data Frame for 1D Covariates
#'
#' Filters a long-format data frame based on LASSO regression, removing cell_ids based on a quasipoisson model. 
#'
#' @param long_df A data frame in long format. The data frame should contain at least 
#'   the columns \code{treatment}, \code{trtlr_cell_id}, and \code{value}. The function
#'   will filter for rows where \code{treatment} is either "ctrl" or "trt", and \code{value}
#'   falls between 10 and a specified cap value (available in the global environment as \code{cap_val}).
#' @param mc.cores The number of cores to use for parallel computation during cross-validation in LASSO.
#'   Default is 1.
#' @param return_df Logical; if \code{TRUE}, the function returns a data frame of retained
#'   rows. If \code{FALSE} (default), it returns a vector of \code{trtlr_cell_id}s that
#'   were removed by the filtering process.
#'
#' @details
#' The function first filters the input data frame based on specific conditions: it retains
#' only rows where \code{treatment} is "ctrl" or "trt", and where \code{value} is greater than 10
#' and less than a global \code{cap_val}. It renames the column \code{trtlr_cell_id} to \code{cell_id} 
#' for the LASSO analysis. A LASSO regression with cross-validation is then applied using the 
#' quasipoisson family. The retained covariates (non-zero beta coefficients) are identified, 
#' and either the filtered data frame or the list of removed covariates is returned, depending on 
#' the \code{return_df} parameter.
#'
#' @return If \code{return_df = TRUE}, a filtered version of the input data frame is returned.
#'   If \code{return_df = FALSE} (default), a vector of removed \code{trtlr_cell_id}s is returned.
#'
#' 
#' @examples
#' \dontrun{
#'   filtered_df <- lasso_1d_filter(my_long_df, mc.cores = 4, return_df = TRUE)
#'   removed_ids <- lasso_1d_filter(my_long_df, mc.cores = 2, return_df = FALSE)
#' }
#' 
#' @export
lasso_1d_filter=function(long_df,mc.cores=1,return_df=F,thresh=-Inf,cap=Inf){
  #browser()
  long_df %>% dplyr::filter(treatment %in% c("ctrl","trt")) %>%  dplyr::filter(value>thresh,value<cap)  %>% dplyr::rename(cell_id=trtlr_cell_id) %>% dplyr::select(treatment,cell_id,value) ->long_df_slim
  X <- model.matrix(value ~ cell_id, data = long_df_slim)[, -1]
  y=long_df_slim$value
  doMC::registerDoMC(cores=mc.cores)
  glmnet::cv.glmnet(X,y,trace.it = T,nfolds=10,gamma=1,parallel=T,family=quasipoisson())->cvfit
  glmnet::glmnet(x = X,y = y,alpha = 1,family=quasipoisson,lambda = cvfit$lambda.1se)->qfit
  qfit$beta %>% as.matrix() %>% as.data.frame() %>% tibble::rownames_to_column(var="cell_id")  %>% dplyr::rename(beta=s0) %>% dplyr::filter(beta!=0,cell_id!="(Intercept)")  %>% dplyr::pull(cell_id) ->retained_cell_ids
  #browser()
  setdiff(long_df$trtlr_cell_id,gsub(x = retained_cell_ids,pattern="cell_id",replacement="")) ->removed_cell_ids
  if(return_df){
    return(long_df %>% dplyr::filter(trtlr_cell_id %in% retained_cell_ids))
  }
  return(removed_cell_ids)
}