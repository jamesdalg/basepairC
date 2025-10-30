#' Smooth and Downsample Hi-C Matrices with MNase Bias Correction
#'
#' This function smooths and downsamples control and treatment matrices along with a MNase bias matrix.
#'
#' @param ctrl_matrix A numeric matrix representing the control interaction counts.
#' @param trt_matrix A numeric matrix representing the treatment interaction counts.
#' @param mnase_bias_matrix A numeric matrix representing the MNase-seq bias, corresponding to the control and treatment matrices.
#' @param smooth_window An integer specifying the window size for smoothing the MNase bias matrix using a kernel smoothing function. Default is `100`.
#' @param ds_factor An integer specifying the downsampling factor to reduce the resolution of the matrices. Default is `10`.
#' @param filter_action A character vector specifying the type of filter to apply to the MNase bias matrix. Options include `mean` (default, using a kernel smooth) and `max` (using dilation).
#' @return A list containing three matrices:
#' \item{ctrl}{The smoothed and downsampled control matrix.}
#' \item{trt}{The smoothed and downsampled treatment matrix.}
#' \item{mnase_bias}{The smoothed and downsampled MNase bias matrix.}
#'
#' @details 
#' The function first sets the row and column names of the MNase bias matrix to match those of the control matrix. It then performs a kernel smoothing operation on the MNase bias matrix to account for proximal effects of MNase-seq bias. After smoothing, the function downsamples the control, treatment, and MNase bias matrices using the specified downsampling factor.
#'
#' @examples
#' \dontrun{
#' smoothed_downsampled_data <- ds_smooth_mat_set(ctrl_matrix, trt_matrix, mnase_bias_matrix)
#' ctrl_matrix_smoothed <- smoothed_downsampled_data$ctrl
#' trt_matrix_smoothed <- smoothed_downsampled_data$trt
#' mnase_bias_smoothed <- smoothed_downsampled_data$mnase_bias
#' }
#'
#' @importFrom smoothie kernel2dsmooth
#' @importFrom dplyr mutate
#' @export
ds_smooth_mat_set<-function(ctrl_matrix,trt_matrix,mnase_bias_matrix,smooth_window=1,ds_factor=10,filter_action=c("mean","max") ,mean_ds=F){
  message_parallel("setting column and row names")
  #add proper row and columna names to the mnase bias matrix
  colnames(mnase_bias_matrix) <- colnames(ctrl_matrix)
  rownames(mnase_bias_matrix) <- rownames(ctrl_matrix)
  #smooth the mnase bias matrix (MNase effects are proximal/surrounding bias high points, and high bias points can be directly in the middle of affected signals but not at the same base pair).
  message_parallel("performing kernel smooth on bias matrix")
  #browser()
  if(filter_action=="mean"){
  smoothie::kernel2dsmooth(mnase_bias_matrix,n=smooth_window,kernel.type="boxcar")->mnase_bias_matrix_smoothed
  }else if(filter_action=="max"){
    EBImage::dilate(mnase_bias_matrix, EBImage::makeBrush(smooth_window, shape="box")) ->mnase_bias_matrix_smoothed
  }
  colnames(mnase_bias_matrix_smoothed) <- colnames(ctrl_matrix)
  rownames(mnase_bias_matrix_smoothed) <- rownames(ctrl_matrix)
  message_parallel("downsampling bias matrix")
  mnase_bias_matrix_smoothed  %>% downsample_with_names(factor=ds_factor,mean=mean_ds) -> mnase_bias_matrix_smoothed_ds
  message_parallel("downsampling ctrl matrix")
  downsample_with_names(ctrl_mat_no_scn_deseq2_diff_ds,factor = ds_factor,mean = mean_ds)->ctrl_ds_mat
  message_parallel("downsampling trt matrix")
  downsample_with_names(trt_mat_no_scn_deseq2_diff_ds,factor = ds_factor,mean=mean_ds)->trt_ds_mat
  message_parallel("finished downsampling!")
  return(list(ctrl=ctrl_ds_mat,trt=trt_ds_mat,mnase_bias=mnase_bias_matrix_smoothed_ds))
}
