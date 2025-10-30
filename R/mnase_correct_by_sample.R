#' Correct MNase Bias by Sample Using Generalized Additive Models (GAM)
#'
#' This function applies a correction for MNase bias in a dataset by fitting a GAM model 
#' for each sample separately and adjusting the values accordingly. The correction is 
#' based on the `mnase_bias_6bp` and `mnase_bias_8bp` covariates in the data.
#'
#' @param cell_df A data.frame or data.table containing the dataset to be corrected. 
#' The dataset must include at least the columns `sample`, `value`, `mnase_bias_6bp`, and `mnase_bias_8bp`.
#' @param single_model_cores An integer specifying the number of CPU cores to use for fitting each model.
#' Defaults to 1 (single-threaded).
#'
#' @return A data.table with the same structure as `cell_df`, but with an additional 
#' column `snr_filtered_value_gam`, which contains the bias-corrected values.
#'
#' @details 
#' This function splits the input dataset by sample, fits a GAM model to the values in 
#' each sample based on the `mnase_bias_6bp` and `mnase_bias_8bp` variables, and computes 
#' a signal-to-noise ratio (SNR) filtered value for each observation. The corrected value is 
#' computed as the ratio of the original value to the predicted value from the GAM model. If 
#' the ratio is greater than 1, the original value is retained; otherwise, the value is set to 0.
#'
#' The `mgcv::gam()` function is used to fit the model, and a negative binomial family (`mgcv::nb()`) 
#' is assumed. The correction is done on a per-sample basis, and the computation is parallelized 
#' using the `pbapply::pblapply()` function.
#'
#' @import mgcv
#' @import data.table
#' @import pbapply
#' @importFrom stats predict
#' 
#' @examples
#' \dontrun{
#'   corrected_df <- mnase_correct_by_sample(cell_df, single_model_cores = 4)
#' }
#'
#' @export
mnase_correct_by_sample=function(cell_df,single_model_cores=1){
  samples=unique(cell_df$sample)
  #
  # if(length(samples)==1){
  #   browser()
  # }
  pbapply::pblapply(1:length(samples),function(i){
    cell_df_sample=cell_df[sample == samples[i]]
    system.time({gam_model_sample <- mgcv::gam(value ~   mnase_bias_6bp +mnase_bias_8bp, #+ te(start1,start2,bs="cs")
                                               data = cell_df_sample,
                                               family = mgcv::nb(),control = list(nthread=single_model_cores))})
    cell_df_sample[, snr_filtered_value_gam:=ifelse(cell_df_sample$value/predict(gam_model_sample,type="response")>1,cell_df_sample$value/predict(gam_model_sample,type="response"),0)]
    return(cell_df_sample)
  }) %>% rbindlist()->cell_df_sample_complete
  return(cell_df_sample_complete)
}