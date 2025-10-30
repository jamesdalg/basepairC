#' @title Select Optimal Lambda Value
#' @description Selects the optimal lambda value for the fused Lasso regression using the fusedLattice method. It computes various statistics and generates plots to aid in selecting the lambda that best fits the data.
#' @param trt_ratio_mat A numeric matrix representing the treatment ratio matrix.
#' @param lambda A numeric vector of lambda values to test. If NULL (default), the function will use internally generated lambda values.
#' @param eps A numeric value specifying the convergence threshold for the fusedLattice function. Default is \code{1e-4}.
#' @param maxIter An integer specifying the maximum number of iterations for the fusedLattice function. Default is \code{1e4}.
#' @param thresh A numeric threshold value for filtering lasso values in plotting. Default is \code{0.1}.
#' @return A numeric value representing the selected optimal lambda.
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Applies the fusedLattice function over a range of lambda values on the provided treatment ratio matrix.
#'   \item Calculates the number of non-zero coefficients in the beta estimates for each lambda.
#'   \item Computes lasso test values and summarizes them using statistical measures such as mean, median, standard deviation, kurtosis, and skewness.
#'   \item Generates various plots including density ridges, histograms, and line plots of statistical summaries against lambda.
#'   \item Determines inflection points in the mean lasso test values versus lambda curve using the \code{inflection} package.
#'   \item Selects the optimal lambda based on the inflection points, typically using the median of the lambda values corresponding to these points.
#' }
#' @examples
#' \dontrun{
#' # Example usage:
#' library(glmgen)
#' library(pbapply)
#' library(dplyr)
#' library(ggplot2)
#' library(ggridges)
#' library(inflection)
#' library(e1071)
#' library(ggrepel)
#'
#' # Create a treatment ratio matrix
#' trt_ratio_mat <- matrix(rnorm(100), nrow = 10)
#'
#' # Select the optimal lambda
#' optimal_lambda <- select_lambda(trt_ratio_mat)
#' }
#' @importFrom glmgen fusedLattice
#' @importFrom pbapply pblapply
#' @importFrom dplyr %>% group_by summarise mutate filter pull
#' @importFrom ggplot2 ggplot aes geom_point geom_line geom_abline labs
#' @importFrom ggridges stat_density_ridges theme_ridges
#' @importFrom inflection ese check_curve bese
#' @importFrom e1071 kurtosis skewness
#' @importFrom ggrepel geom_text_repel
#' @export
select_lambda=function(trt_ratio_mat,lambda=NULL,eps=1e-4,maxIter=1e4,thresh=0.1){  
  #lambda=seq(-5,40,1)
  glmgen::fusedLattice(as.matrix(trt_ratio_mat), lambda = lambda, eps=eps, maxIter=maxIter)->lambda_test_res
  #for every value of lambda, find the number of nonzero points in beta.
  pbapply::pblapply(1:length(lambda_test_res$lambda),function(i) {
    #sum(glmgen::fusedLattice(as.matrix(trt_ratio_mat), lambda = lambda_test_res$lambda[i], eps=eps, maxIter=maxIter)$beta!=0)
    if(i==1){return(
      length(lambda_test_res$beta[1:(nrow(trt_ratio_mat)*ncol(trt_ratio_mat)*(i) )]!=0)
    )}
    if(i>1){return(length(
      lambda_test_res$beta[(nrow(trt_ratio_mat)*nrow(trt_ratio_mat)*(i-1)):(nrow(trt_ratio_mat)*ncol(trt_ratio_mat)*(i) )]
      !=0))}
  }) %>% unlist()->lambda_test_res_nonzero
  
  pbapply::pblapply(1:length(lambda_test_res$lambda),function(i) {
    #foreach(i=1:length(lambda_test_res$lambda),.combine="c") %do% {
    #sum(glmgen::fusedLattice(as.matrix(trt_ratio_mat), lambda = lambda_test_res$lambda[i], eps=eps, maxIter=maxIter)$beta!=0)
    if(i==1){
      #browser()
      return(
        
        as.numeric(matrix((lambda_test_res$beta[1:(nrow(trt_ratio_mat)*ncol(trt_ratio_mat)*(i))]),nrow=nrow(trt_ratio_mat))*trt_ratio_mat)
      )}
    if(i>1){
      #browser()
      #print(nrow(trt_ratio_mat) * nrow(trt_ratio_mat) * (i-1) )
      #print(paste0(((nrow(trt_ratio_mat)*nrow(trt_ratio_mat)*(i-1))-1))," to ",(nrow(trt_ratio_mat)*ncol(trt_ratio_mat)*(i) ))
      return(
        
        as.numeric(trt_ratio_mat*matrix(lambda_test_res$beta[((nrow(trt_ratio_mat)*nrow(trt_ratio_mat)*(i-1))):(nrow(trt_ratio_mat)*ncol(trt_ratio_mat)*(i) -1)],nrow=nrow(trt_ratio_mat)))
      )}
  }) ->lambda_lasso_test_vals
  #put each set into a data frame with the lambda value and all the values of the lasso test, then cbind them all and plot them as histograms in ggplot, using a ggridges geom.
  lapply(1:length(lambda_test_res$lambda),function(i) {
    return(
      data.frame(lambda=lambda_test_res$lambda[i],lasso_vals=lambda_lasso_test_vals[[i]])
    )
  }) %>% do.call(rbind,.)->lambda_lasso_test_df
  library(ggplot2)
  library(ggridges)
  
  #browser()
  ggplot(lambda_lasso_test_df,aes(x=lasso_vals,y=lambda,fill=lambda,group=lambda)) + stat_density_ridges(geom="density_ridges_gradient",calc_ecdf=T) + theme_ridges() + scale_fill_viridis_c() + labs(title="Lasso Test Values for Different Lambda Values",x="Lasso Test Values",y="Lambda Value")
  lambda_lasso_test_df %>% dplyr::group_by(lambda) %>% dplyr::summarise(mean_lasso=mean(lasso_vals),median_lasso=median(lasso_vals),sd_lasso=sd(lasso_vals),first_quartile=quantile(lasso_vals,.25),third_quartile=quantile(lasso_vals,.75),n_pos=sum(which(lasso_vals>0)),kurtosis=e1071::kurtosis(lasso_vals),skew=e1071::skewness(lasso_vals))->lambda_lasso_test_df_summ
  plot(lambda_lasso_test_df_summ$mean_lasso,x=lambda_lasso_test_df_summ$lambda,type="b",xlab="Lambda",ylab="Mean Lasso Test Value",main="Mean Lasso Test Value vs Lambda")
  #do kurtosis
  plot(lambda_lasso_test_df_summ$kurtosis,x=lambda_lasso_test_df_summ$lambda,type="b",xlab="Lambda",ylab="Kurtosis Lasso Test Value",main="Kurtosis Lasso Test Value vs Lambda")
  #do skewness
  plot(lambda_lasso_test_df_summ$skew,x=lambda_lasso_test_df_summ$lambda,type="b",xlab="Lambda",ylab="Skewness Lasso Test Value",main="Skewness Lasso Test Value vs Lambda")
  inflection::ese(lambda_lasso_test_df_summ$mean_lasso,lambda_lasso_test_df_summ$lambda,inflection::check_curve(lambda_lasso_test_df_summ$mean_lasso,lambda_lasso_test_df_summ$lambda)$index)[1:2]->inflection_points
  inflection_points=unique(c(which(diff(sign(diff(diff(lambda_lasso_test_df_summ$mean_lasso))))!=0),inflection_points))
  inflection::bese(lambda_lasso_test_df_summ$mean_lasso,lambda_lasso_test_df_summ$lambda,inflection::check_curve(lambda_lasso_test_df_summ$mean_lasso,lambda_lasso_test_df_summ$lambda)$index)$ipblast->ipblast_val
  lambda_lasso_test_df_summ %>% dplyr::mutate(index=1:nrow(.)) %>%
    ggplot(aes(y=mean_lasso,x=lambda,label=ifelse(index %in% inflection_points,paste0(lambda,",IP"),paste0(lambda,",",as.character(index))))) + geom_point() + geom_line() + ggrepel::geom_text_repel() + geom_abline(aes(intercept=ipblast_val),slope=0,linetype="dashed") + labs(title="Mean Lasso Test Values for Different Lambda Values",x="Mean Lasso Test Values",y="Lambda Value") 
  #ggplot version of the same with l
  #plot with standard geom histogram:
  ggplot(lambda_lasso_test_df %>% dplyr::filter(lasso_vals>thresh),aes(x=lasso_vals,fill=lambda)) + geom_histogram(binwidth = 0.1) + facet_wrap(~lambda,scales="free") + labs(title="Lasso Test Values for Different Lambda Values",x="Lasso Test Values",y="Frequency") + geom_density() 
  #browser()
  #plot(lambda_test_res$lambda,main="fusedLattice generated Lambda values",ylab="Lambda",ylab="Index")
  lambda=lambda_lasso_test_df_summ %>% dplyr::mutate(index=1:nrow(.)) %>% dplyr::filter(index %in% inflection_points) %>% dplyr::pull(lambda) %>% median() #try median if this doens't work well enough.
  return(lambda)
}