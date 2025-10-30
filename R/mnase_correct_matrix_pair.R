#' Correct MNase Bias in Hi-C Matrices for Control and Treatment
#'
#' This function corrects the MNase-seq bias in two Hi-C matrices (control and treatment) using a pre-computed MNase bias matrix.
#'
#' @param ctrl_mat_ds A numeric matrix representing the control Hi-C interaction counts.
#' @param trt_mat_ds A numeric matrix representing the treatment Hi-C interaction counts.
#' @param mnase_bias_matrix_sm_ds A numeric matrix representing the MNase-seq bias. This matrix should have the same dimensions as `ctrl_mat_ds` and `trt_mat_ds`.
#' @param snr_filter A logical value indicating whether to apply signal-to-noise ratio (SNR) filtering. Default is `FALSE`.
#' @param two_stage A logical value indicating whether to use a two-stage correction. Default is `FALSE`.
#' @param single_model_cores An integer specifying the number of cores to use for the single model correction. Default is `1`.
#' @param residual Does a correction with the residuals instead of outcome. Stronger SNR filter, though the residuals may have a longer tailed distribution. Default is `FALSE`.
#' @return A list containing two corrected Hi-C matrices: 
#' \item{ctrl}{Corrected control Hi-C matrix.}
#' \item{trt}{Corrected treatment Hi-C matrix.}
#'
#' @details 
#' The function first melts the MNase bias matrix and joins it with the melted control and treatment matrices. It then applies a Negative Binomial generalized linear model (GLM) to correct for MNase bias. The correction factor is calculated based on the GLM coefficients and applied to the original interaction counts to produce the corrected matrices. The function outputs the corrected matrices in a list.
#'
#' @examples
#' \dontrun{
#' corrected_matrices <- mnase_correct_matrix_pair(ctrl_mat_ds, trt_mat_ds, mnase_bias_matrix_sm_ds)
#' corrected_ctrl <- corrected_matrices$ctrl
#' corrected_trt <- corrected_matrices$trt
#' }
#'
#' @importFrom dplyr inner_join mutate rename group_by summarize ungroup select filter distinct arrange
#' @importFrom tidyr pivot_longer
#' @importFrom reshape2 melt
#' @importFrom MASS glm.nb
#' @importFrom tibble tibble
#' @export
mnase_correct_matrix_pair<-function(ctrl_mat_ds,trt_mat_ds,mnase_bias_matrix_sm_ds,snr_filter=F,two_stage=F,single_model_cores=1,residual=F){
  #citation: Jin Li, p226-228, Spatial Predictive Modeling with R, 2022.
  #input: two matices, one for control and one for treatment.
  #output: two matrices, one for control and one for treatment, with the same dimensions.
  mnase_bias_matrix_sm_ds %>% reshape2::melt() ->mnase_bias_matrix_melt
  #test this to ensure that the columns are being named correctly.
  #browser()
  mnase_bias_matrix_melt %>% as.data.frame() %>% dplyr::inner_join(ctrl_mat_ds %>% reshape2::melt(),by=c("Var1","Var2"),suffix = c("_mnase",""),copy=T) %>%
    dplyr::inner_join(trt_mat_ds  %>% reshape2::melt(),by=c("Var1","Var2"),suffix = c("_ctrl","_trt"),copy=T) %>%
    dplyr::rename(loc1=Var1,loc2=Var2,ctrl=value_ctrl,trt=value_trt)  %>%
    tidyr::pivot_longer(cols=c("ctrl","trt"),names_to="treatment",values_to = "value") %>% 
    dplyr::rename(mnase=value_mnase)->merged_ctrl_trt_mnase
  
  total_reads_df=tibble::tibble(treatment=c("ctrl","trt"),total_reads=c(sum(ctrl_mat_ds),sum(trt_mat_ds)))
  print(total_reads_df)

  if(two_stage){ #two stage correction
    #this extends the original one-stage correction and uses a nonlinear gam and linear glm to correct for MNase bias.
    total_reads_df=tibble::tibble(treatment=c("ctrl","trt"),total_reads=c(sum(ctrl_mat_ds),sum(trt_mat_ds))) %>% dplyr::mutate(tr_ratio=1/(total_reads/min(total_reads)))
    message_parallel("step 1: fitting inital adaptive smooth on MNase")
    merged_ctrl_trt_mnase %>% dplyr::inner_join(total_reads_df) %>% dplyr::mutate(mnase=basepairC::ztrans(log(mnase+1)))->merged_ctrl_trt_mnase_trans
    mgcv::bam(value~s(mnase,bs="ad"),data=merged_ctrl_trt_mnase_trans,
        family=mgcv::nb(link="log"),method="fREML",control=list(nthread=single_model_cores,ncv.threads=single_model_cores),nthreads=single_model_cores,discrete=T)->mnase_gam
    message_parallel("step 2: fitting a linear mnase fit, using smoothed model residuals from step 1 as a covariate ") #residuals to NB GLM, fully linear
    

    resid_data_df=merged_ctrl_trt_mnase %>% dplyr::inner_join(total_reads_df) %>% dplyr::mutate(mnase=basepairC::ztrans(log(mnase+1))) %>% dplyr::mutate(fitted=predict(mnase_gam,type="response"),residual=value-fitted)
    #browser()
    #glm.nb(value~residual+mnase,data=resid_data_df)->mnase_glm_resid
    #biglm::bigglm(value~residual+mnase,data=resid_data_df,family=quasipoisson)->mnase_glm_resid
    #THIS WORKS (BELOW):
    mgcv::bam(value~s(residual,bs="ad")+mnase+treatment,data=resid_data_df,family=mgcv::nb(link="log"),method="fREML",control=list(nthread=single_model_cores,mgcv.half=30,maxit=1e3,trace=T),nthreads=single_model_cores,samfrac=.001,use.chol = T,discrete=T)->mnase_glm_resid
    #glmmTMB::glmmTMB(value~residual+mnase,data=resid_data_df,family=poisson(),control=glmmTMB::glmmTMBControl(parallel=parallel::detectCores(),optimizer=optim,collect=F,profile=T,optArgs=list(method="CG")))->mnase_glm_resid_pois
    #THIS ALSO WORKS:
    #glmmTMB::glmmTMB(value~s(residual)+mnase+treatment,data=resid_data_df,family=glmmTMB::nbinom2(),control=glmmTMB::glmmTMBControl(parallel=parallel::detectCores(),optimizer=optim,collect=F,profile=T,optCtrl = list(trace=1)),REML=T)->mnase_glm_resid_nb2 
    #mgcv::bam(value~residual+mnase,data=resid_data_df,family=mgcv::nb(link="log"),method="fREML",control=list(nthread=single_model_cores),nthreads=single_model_cores,discrete=T)->mnase_glm_resid
  if(snr_filter){  #two stage, with snr
    message_parallel("Correcting with SNR filter")
    dplyr::inner_join(merged_ctrl_trt_mnase,total_reads_df,
                      by=c("treatment"))  %>%
    dplyr::mutate(correction_factor=exp(coef(mnase_glm_resid)["(Intercept)"] +
                                      coef(mnase_glm_resid)["mnase"] * mnase ),
    corrected=ifelse(value/correction_factor>1,value/correction_factor,0)*tr_ratio
    ) ->long_df_coarse_w_results_corrected    
  } else { #two stage, no snr
    message_parallel("Correcting without SNR filter")
    dplyr::inner_join(merged_ctrl_trt_mnase,total_reads_df,
                      by=c("treatment"))  %>%
      dplyr::mutate(correction_factor=exp(coef(mnase_glm_resid)["(Intercept)"] +
                                            coef(mnase_glm_resid)["mnase"] * mnase ),
                    corrected=(value/correction_factor)*tr_ratio
      ) ->long_df_coarse_w_results_corrected
  }  
    
  } else{ #single stage correction
    #this follows HiCNorm's approach, but uses the NB distribution, which is more appropriate for skewed count data.
    #browser()
    merged_data=merged_ctrl_trt_mnase %>% dplyr::inner_join(total_reads_df) %>% dplyr::mutate(min_mnase=min(nonzero(mnase)),mnase=basepairC::ztrans(log(mnase+min_mnase)))
  MASS::glm.nb(round(value)~mnase+treatment,data=merged_data ) ->mnase_glm  #+treatment #+offset(log(total_reads)
  print(summary(mnase_glm))
  print("mnase*beta_mnase distribution:")
  print((coef(mnase_glm)["mnase"]*(merged_ctrl_trt_mnase %>% dplyr::inner_join(total_reads_df) %>% dplyr::mutate(min_mnase=min(mnase),mnase=basepairC::ztrans(log(mnase+min_mnase))) %>% dplyr::pull(mnase))  ) %>% summary())
  #browser()
  if(residual){
    if(snr_filter){
      #simplify notation: (min(total_reads)/total_reads) is a simpler way of writing things.
      dplyr::inner_join(merged_ctrl_trt_mnase,merged_ctrl_trt_mnase %>% dplyr::group_by(treatment) %>% dplyr::summarize(total_reads=sum(value)) %>% dplyr::ungroup() %>% dplyr::mutate(tr_ratio=1/(total_reads/min(total_reads))),by=c("treatment"))   %>%   dplyr::mutate(correction_factor=exp(coef(mnase_glm)[1] + coef(mnase_glm)[2] * mnase ),
                                                                                                                                                                                                                                                                          corrected=ifelse(residuals(mnase_glm)/correction_factor>1,residuals(mnase_glm),0)*tr_ratio
      ) ->long_df_coarse_w_results_corrected} else{
        dplyr::inner_join(merged_ctrl_trt_mnase,merged_ctrl_trt_mnase %>% dplyr::group_by(treatment) %>% dplyr::summarize(total_reads=sum(value)) %>% dplyr::ungroup() %>% dplyr::mutate(tr_ratio=1/(total_reads/min(total_reads))),by=c("treatment"))  %>%   dplyr::mutate(correction_factor=exp(coef(mnase_glm)[1] + coef(mnase_glm)[2] * mnase ),
                                                                                                                                                                                                                                                                            corrected=(residuals(mnase_glm))*tr_ratio
        ) ->long_df_coarse_w_results_corrected
      }
  
  
  }
  if(!residual){
  if(snr_filter){
    #simplify notation: (min(total_reads)/total_reads) is a simpler way of writing things.
  dplyr::inner_join(merged_ctrl_trt_mnase,merged_ctrl_trt_mnase %>% dplyr::group_by(treatment) %>% dplyr::summarize(total_reads=sum(value)) %>% dplyr::ungroup() %>% dplyr::mutate(tr_ratio=1/(total_reads/min(total_reads))),by=c("treatment"))  %>%   dplyr::mutate(correction_factor=exp(coef(mnase_glm)[1] + coef(mnase_glm)[2] * mnase ),
                                                                                                                                                                                                                                                                      corrected=ifelse(value/correction_factor>1,value/correction_factor,0)*tr_ratio
  ) ->long_df_coarse_w_results_corrected} else{
    dplyr::inner_join(merged_ctrl_trt_mnase,merged_ctrl_trt_mnase %>% dplyr::group_by(treatment) %>% dplyr::summarize(total_reads=sum(value)) %>% dplyr::ungroup() %>% dplyr::mutate(tr_ratio=1/(total_reads/min(total_reads))),by=c("treatment"))  %>%   dplyr::mutate(correction_factor=exp(coef(mnase_glm)[1] + coef(mnase_glm)[2] * mnase ),
                                                                                                                                                                                                                                                                        corrected=(value/correction_factor)*tr_ratio
    ) ->long_df_coarse_w_results_corrected
  }
  }
  }
  message_parallel("finalizing results")
  #same both both one stage and two stage
  long_df_coarse_w_results_corrected %>% dplyr::select(loc1,loc2,treatment,corrected) %>% dplyr::filter(treatment=="ctrl") %>% dplyr::select(loc1,loc2,corrected) %>% dplyr::mutate(loc1=as.character(loc1),loc2=as.character(loc2))  %>% dplyr::distinct() %>% dplyr::arrange(loc1,loc2)   %>%  freeze(loc1_col = "loc1",loc2_col="loc2",value_col = "corrected")->ctrl_ds_corrected
  long_df_coarse_w_results_corrected %>% dplyr::select(loc1,loc2,treatment,corrected) %>% dplyr::filter(treatment=="trt") %>% dplyr::select(loc1,loc2,corrected) %>% dplyr::mutate(loc1=as.character(loc1),loc2=as.character(loc2))  %>% dplyr::distinct() %>% dplyr::arrange(loc1,loc2)  %>%  freeze(loc1_col = "loc1",loc2_col="loc2",value_col = "corrected")->trt_ds_corrected
  return(list(ctrl=ctrl_ds_corrected,trt=trt_ds_corrected,mnase_glm=mnase_glm,model_data=merged_data, merged_ctrl_trt_mnase=merged_ctrl_trt_mnase,mnase_glm=mnase_glm
              ))
  
}