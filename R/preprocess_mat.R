#' @title preprocess_mat
#' @description This function preprocesses a matrix by applying fast combining of , outlier removal, and fused 2D lasso filtering.
#' @param fn The file name of the matrix to be preprocessed.
#' @param cond The condition of the matrix.
#' @param bin_size The bin size of the matrix.
#' @param thresh Removes reads below this count value from the input matrix
#' @param cap Reads above this height will be discarded or thresholded.
#' @param lambda The lambda value for the fused 2D lasso.
#' @param eps The epsilon value for the fused 2D lasso.
#' @param maxIter The maximum number of iterations for the fused 2D lasso.
#' @param quantile The quantile value for the cap.
#' @param zero_top Logical, indicating whether to zero out top values (if 0, then the cap value is used).
#' @param chr The chromosome of the matrix.
#' @param start_row The starting row of the matrix.
#' @param end_row The ending row of the matrix.
#' @param start_col The starting column of the matrix.
#' @param end_col The ending column of the matrix.
#' @param dsfactor The downsample factor for the matrix (essentially a 6000x6000 matrix becomes a 2000x2000 matrix if the dsfactor is 3).
#' @param mat_start_row The starting row of the matrix.
#' @param mat_end_row The ending row of the matrix.
#' @param mat_start_col The starting column of the matrix.
#' @param mat_end_col The ending column of the matrix.
#' @param plot Logical, indicating whether to plot the matrix.
#' @param orig_seg_max_col The maximum column value for the original segmentation.
#' @param orig_seg_max_row The maximum row value for the original segmentation.
#' @param lasso_seg_max_row The maximum row value for the lasso segmentation.
#' @param lasso_seg_max_col The maximum column value for the lasso segmentation.
#' @param orig_seg_method The distribution for the original segmentation (P, G, or B for poisson, gaussian, or binomial, respectively).
#' @param lasso_seg_method The distribution for the lasso segmentation.
#' @param orig_seg_model The model for the original matrix segmentation (D or Dplus).
#' @param lasso_seg_model The model for the lasso matrix segmentation.
#' @param debug Logical, indicating whether to print debug information.
#' @param n_cores The number of cores to use for parallel processing.
#' @param make_symmetric Logical, indicating whether to make the matrix symmetric.
#' @param ds_fn The function to use for downsampling.
#' @param seg_algorithm The algorithm to use for segmentation.
#' @return A list containing the processed and unprocessed matrices.
#' @importFrom foreach foreach `%dopar%` `%do%`
#' @importFrom magrittr `%>%` 
#' @importFrom GenomicRanges GRanges
#' @importFrom Matrix forceSymmetric
#' @export preprocess_mat
preprocess_mat=function(fn,cond,bin_size=1,thresh=0,cap=Inf,lambda=5,eps=1e-10,maxIter=20,quantile=0,zero_top=F,chr="chr15",start_row=61982445,end_row=61988445,start_col=61982445,end_col=61988445,dsfactor=50,
                        mat_start_row=1,mat_end_row=6000,mat_start_col=1,mat_end_col=6000,plot=F,orig_seg_max_col=9,orig_seg_max_row=9,lasso_seg_max_row=9,lasso_seg_max_col=9,
                        orig_seg_method="G",lasso_seg_method="G",orig_seg_model="Dplus",lasso_seg_model="Dplus",debug=F,n_cores=1,make_symmetric=T,ds_fn="sum",seg_algorithm="jointSeg",
                        hicseg_bp_selection=T,lasso2d=F,verbose=F,segments_orig_obj=NULL,segments_lasso_obj=NULL,seg_method="RBS",smooth=T,smooth_window=5,lasso_seg=F) {
  thresh=thresh*dsfactor
  options(nThread=n_cores)
  if(debug){browser()}
  if(verbose){message_parallel('reading data')}
  if(any(class(fn) %in% "matrix")){
    if(make_symmetric){
      orig_mat_norm=(fn %>% basepairC:::make_symmetric() %>% as("sparseMatrix")) %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col]}else{
        orig_mat_norm=(fn %>% as.matrix()  %>% as("sparseMatrix")) %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col]  
      }  
  }
  if(any(class(fn) %in% "dgCMatrix")){
    if(make_symmetric){
      orig_mat_norm=fn %>% basepairC:::make_symmetric()  %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col]}else{
        orig_mat_norm=fn %>% as.matrix()   %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col]
      }
  }
  if(any(class(fn) %in% "character")){
    if(make_symmetric){
      orig_mat_norm=(data.table::fread(fn,nThread=n_cores,showProgress = T) %>% as.matrix() %>% basepairC:::make_symmetric() %>% as("sparseMatrix")) %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col]}else{
        orig_mat_norm=(data.table::fread(fn,nThread=n_cores,showProgress = T) %>% as.matrix()  %>% as("sparseMatrix")) %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col]  
      }
  }
  GRanges(seqnames=chr, IRanges(start=seq(from=start_row,to=end_row,by=bin_size), end=seq(from=start_row+bin_size-1,to=end_row+bin_size-1,by=bin_size)))->row_gr
  GRanges(seqnames=chr, IRanges(start=seq(from=start_col,to=end_col,by=bin_size), end=seq(from=start_col+bin_size-1,to=end_col+bin_size-1,by=bin_size)))->col_gr
  underscored_positions_col=CNVScope::GRanges_to_underscored_pos(col_gr)
  underscored_positions_row=CNVScope::GRanges_to_underscored_pos(row_gr)
  col_gr %>% CNVScope::GRanges_to_underscored_pos() %>% .[mat_start_col:(mat_end_col)] ->colnames(orig_mat_norm)
  row_gr %>% CNVScope::GRanges_to_underscored_pos() %>% .[mat_start_row:(mat_end_row)]->rownames(orig_mat_norm)
  #browser()
  #n_cores=parallel::detectCores()/2
  if(verbose){message_parallel("downsampling")}
  if(dsfactor!=1){
  CNVScope::downsample_genomic_matrix(as.matrix(orig_mat_norm),downsamplefactor = dsfactor,singlechromosome = T)->orig_mat_norm_small
  if(ds_fn=="sum"){orig_mat_norm_small_pooled <- block_sum_matrix(orig_mat_norm, dsfactor, dsfactor) %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)} 
  if(ds_fn=="mean"){orig_mat_norm_small_pooled <- block_sum_matrix(orig_mat_norm, dsfactor, dsfactor,average = T) %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)}
  colnames(orig_mat_norm_small)->colnames(orig_mat_norm_small_pooled)
  rownames(orig_mat_norm_small)->rownames(orig_mat_norm_small_pooled)} else{
    orig_mat_norm_small=as.matrix(orig_mat_norm)
    orig_mat_norm_small_pooled=as.matrix(orig_mat_norm)
  }
  if(verbose){message_parallel("2D lasso cleaning")}
 # if(lasso2d){lambda=1}
  if(lasso2d){
  # if(is.null(lambda)){
  #   glmgen::fusedLattice(as.matrix(orig_mat_norm_small_pooled), lambda = lambda, eps=eps, maxIter=maxIter)->lambda_test_res
  #   #for every value of lambda, find the number of nonzero points in beta.
  #   pbapply::pblapply(1:length(lambda_test_res$lambda),function(i) {
  #     if(i==1){return(
  #       length(lambda_test_res$beta[1:(nrow(orig_mat_norm_small_pooled)*ncol(orig_mat_norm_small_pooled)*(i) )]!=0)
  #     )}
  #     if(i>1){return(length(
  #       lambda_test_res$beta[(nrow(orig_mat_norm_small_pooled)*nrow(orig_mat_norm_small_pooled)*(i-1)):(nrow(orig_mat_norm_small_pooled)*ncol(orig_mat_norm_small_pooled)*(i) )]
  #       !=0))}
  #   }) %>% unlist()->lambda_test_res_nonzero
  #   
  #   pbapply::pblapply(1:length(lambda_test_res$lambda),function(i) {
  #     if(i==1){
  #       return(
  #         as.numeric(matrix((lambda_test_res$beta[1:(nrow(orig_mat_norm_small_pooled)*ncol(orig_mat_norm_small_pooled)*(i))]),nrow=nrow(orig_mat_norm_small_pooled))*orig_mat_norm_small_pooled)
  #       )}
  #     if(i>1){
  #       
  #       return(
  #         
  #         as.numeric(orig_mat_norm_small_pooled*matrix(orig_mat_norm_small_pooled$beta[((nrow(orig_mat_norm_small_pooled)*nrow(orig_mat_norm_small_pooled)*(i-1))):(nrow(orig_mat_norm_small_pooled)*ncol(orig_mat_norm_small_pooled)*(i) -1)],nrow=nrow(orig_mat_norm_small_pooled)))
  #       )}
  #   }) ->lambda_lasso_test_vals
  #   #put each set into a data frame with the lambda value and all the values of the lasso test, then cbind them all and plot them as histograms in ggplot, using a ggridges geom.
  #   lapply(1:length(lambda_test_res$lambda),function(i) {
  #     return(
  #       data.frame(lambda=lambda_test_res$lambda[i],lasso_vals=lambda_lasso_test_vals[[i]])
  #     )
  #   }) %>% do.call(rbind,.)->lambda_lasso_test_df
  #   library(ggplot2)
  #   library(ggridges)
  #   ggplot(lambda_lasso_test_df,aes(x=lasso_vals,y=lambda,fill=lambda,group=lambda)) + stat_density_ridges(geom="density_ridges_gradient",calc_ecdf=T) + theme_ridges() + scale_fill_viridis_c() + labs(title="Lasso Test Values for Different Lambda Values",x="Lasso Test Values",y="Lambda Value")
  #   lambda_lasso_test_df %>% dplyr::group_by(lambda) %>% dplyr::summarise(mean_lasso=mean(lasso_vals),median_lasso=median(lasso_vals),sd_lasso=sd(lasso_vals),first_quartile=quantile(lasso_vals,.25),third_quartile=quantile(lasso_vals,.75),n_pos=sum(which(lasso_vals>0)),kurtosis=e1071::kurtosis(lasso_vals),skew=e1071::skewness(lasso_vals))->lambda_lasso_test_df_summ
  #   plot(lambda_lasso_test_df_summ$mean_lasso,x=lambda_lasso_test_df_summ$lambda,type="b",xlab="Lambda",ylab="Mean Lasso Test Value",main="Mean Lasso Test Value vs Lambda")
  #   #do kurtosis
  #   plot(lambda_lasso_test_df_summ$kurtosis,x=lambda_lasso_test_df_summ$lambda,type="b",xlab="Lambda",ylab="Kurtosis Lasso Test Value",main="Kurtosis Lasso Test Value vs Lambda")
  #   #do skewness
  #   plot(lambda_lasso_test_df_summ$skew,x=lambda_lasso_test_df_summ$lambda,type="b",xlab="Lambda",ylab="Skewness Lasso Test Value",main="Skewness Lasso Test Value vs Lambda")
  #   inflection::ese(lambda_lasso_test_df_summ$mean_lasso,lambda_lasso_test_df_summ$lambda,inflection::check_curve(lambda_lasso_test_df_summ$mean_lasso,lambda_lasso_test_df_summ$lambda)$index)[1:2]->inflection_points
  #   inflection_points=unique(c(which(diff(sign(diff(diff(lambda_lasso_test_df_summ$mean_lasso))))!=0),inflection_points))
  #   inflection::bese(lambda_lasso_test_df_summ$mean_lasso,lambda_lasso_test_df_summ$lambda,inflection::check_curve(lambda_lasso_test_df_summ$mean_lasso,lambda_lasso_test_df_summ$lambda)$index)$ipblast->ipblast_val
  #   lambda_lasso_test_df_summ %>% dplyr::mutate(index=1:nrow(.)) %>%
  #     ggplot(aes(y=mean_lasso,x=lambda,label=ifelse(index %in% inflection_points,paste0(lambda,",IP"),paste0(lambda,",",as.character(index))))) + geom_point() + geom_line() + ggrepel::geom_text_repel() + geom_abline(aes(intercept=ipblast_val),slope=0,linetype="dashed") + labs(title="Mean Lasso Test Values for Different Lambda Values",x="Mean Lasso Test Values",y="Lambda Value") 
  #   #ggplot version of the same with l
  #   #plot with standard geom histogram:
  #   ggplot(lambda_lasso_test_df %>% dplyr::filter(lasso_vals>thresh),aes(x=lasso_vals,fill=lambda)) + geom_histogram(binwidth = 0.1) + facet_wrap(~lambda,scales="free") + labs(title="Lasso Test Values for Different Lambda Values",x="Lasso Test Values",y="Frequency") + geom_density() 
  #   #plot(lambda_test_res$lambda,main="fusedLattice generated Lambda values",ylab="Lambda",ylab="Index")
  #   lambda=lambda_lasso_test_df_summ %>% dplyr::mutate(index=1:nrow(.)) %>% dplyr::filter(index %in% inflection_points) %>% dplyr::pull(lambda) %>% median() #try median if this doens't work well enough.
  #   as.matrix(matrix(minmaxnorm(glmgen::fusedLattice(as.matrix(orig_mat_norm_small_pooled), lambda = lambda, eps=eps, maxIter=maxIter)$beta,nrow=nrow(orig_mat_norm_small_pooled) )  )*orig_mat_norm_small_pooled)->lasso_cleaned_mat
  # } else {
    
    as.matrix(matrix(minmaxnorm(glmgen::fusedLattice(as.matrix(orig_mat_norm_small_pooled), lambda = lambda, eps=eps, maxIter=maxIter)$beta),nrow=nrow(orig_mat_norm_small_pooled)   )*orig_mat_norm_small_pooled)->lasso_cleaned_mat
#  } 
  } else {
    lasso_cleaned_mat=orig_mat_norm_small_pooled
  }
  if(verbose){message_parallel("determining optimal segment breakpoints")}
  if(debug){browser()}
  
  if(seg_algorithm=="HiCseg"&hicseg_bp_selection){
    HiCseg::HiCseg_linkC_R(mat_data=orig_mat_norm_small_pooled,nb_change_max = nrow(orig_mat_norm_small_pooled),distrib=orig_seg_method,size_mat=max(nrow(orig_mat_norm_small_pooled),ncol(orig_mat_norm_small_pooled)),model="Dplus")->hicseg_res_row
    HiCseg::HiCseg_linkC_R(mat_data=orig_mat_norm_small_pooled,nb_change_max = nrow(orig_mat_norm_small_pooled),distrib=orig_seg_method,size_mat=max(nrow(orig_mat_norm_small_pooled),ncol(orig_mat_norm_small_pooled)),model="Dplus")->hicseg_res_col
    #do the same for lasso cleaned matrix
    HiCseg::HiCseg_linkC_R(mat_data=lasso_cleaned_mat,nb_change_max = nrow(lasso_cleaned_mat),distrib=orig_seg_method,size_mat=max(nrow(lasso_cleaned_mat),ncol(lasso_cleaned_mat)),model="Dplus")->hicseg_res_row_lasso
    HiCseg::HiCseg_linkC_R(mat_data=lasso_cleaned_mat,nb_change_max = nrow(lasso_cleaned_mat),distrib=orig_seg_method,size_mat=max(nrow(lasso_cleaned_mat),ncol(lasso_cleaned_mat)),model="Dplus")->hicseg_res_col_lasso
    which.max(hicseg_res_row$J)->orig_seg_max_row
    which.max(hicseg_res_col$J)->orig_seg_max_col
    which.max(hicseg_res_row_lasso$J)->lasso_seg_max_row
    which.max(hicseg_res_col_lasso$J)->lasso_seg_max_col
    
    if(verbose){message_parallel("HiCseg optimal segments (row):")}
    if(verbose){message_parallel(orig_seg_max_row)}
    if(verbose){message_parallel("HiCseg optimal segments (col):")}
    if(verbose){message_parallel(orig_seg_max_col)}
    if(verbose){message_parallel("2D LASSO HiCseg optimal segments (row):")}
    if(verbose){message_parallel(lasso_seg_max_row)}
    if(verbose){message_parallel("2D LASSO HiCseg optimal segments (col):")}
    if(verbose){message_parallel(lasso_seg_max_col)}
    #plot hicseg_res_row$J vs hicseg_res_row$t_hat
    if(plot){
      browser()
      ggplot2::ggplot(hicseg_res_row, aes(x = t_hat, y = J)) + ggplot2::geom_line() + ggplot2::geom_vline(xintercept =  hicseg_res_row$t_hat[orig_seg_max_row], color = "red") + ggplot2::labs(title = "HiCseg J vs t_hat (row)", x = "t_hat", y = "J")
      #do the same for the lasso cleaned matrix
      ggplot2::ggplot(hicseg_res_row_lasso, aes(x = t_hat, y = J)) + ggplot2::geom_line() + ggplot2::geom_vline(xintercept =  hicseg_res_row_lasso$t_hat[lasso_seg_max_row], color = "red") + ggplot2::labs(title = "HiCseg J vs t_hat (row) (Lasso Cleaned)", x = "t_hat", y = "J")
    }
  } else{
    if(is.null(orig_seg_max_row)){
      if(nrow(orig_mat_norm_small_pooled)==ncol(orig_mat_norm_small_pooled)){
        #need to test to see if this is much different. If not, definitely just use the J from HiCseg
        orig_near_eigenvalues_row=eigen(Matrix::nearPD(cov(orig_mat_norm_small_pooled))$mat)$values
        orig_near_eigenvalues_col=orig_near_eigenvalues_row
      } else { orig_near_eigenvalues_row=svd(orig_mat_norm_small_pooled)$d;orig_near_eigenvalues_col=svd(gene_ctrl_mat_norm_small_pooled)}
      if(is.null(lasso_seg_max_row)){
        if(nrow(lasso_cleaned_mat)==ncol(lasso_cleaned_mat)){
          lasso_near_eigenvalues_row=eigen(Matrix::nearPD(cov(lasso_cleaned_mat))$mat)$values
          lasso_near_eigenvalues_col=lasso_near_eigenvalues_row
        } else { lasso_near_eigenvalues_row=svd(lasso_cleaned_mat)$d;lasso_near_eigenvalues_col=svd(lasso_cleaned_mat)}
        
        orig_seg_anal_row=nFactors::nScree(orig_near_eigenvalues)
        orig_seg_max_row=floor(median(as.numeric(as.matrix(orig_seg_anal_row$Components))))
        orig_seg_anal_col=nFactors::nScree(orig_near_eigenvalues)
        orig_seg_max_col=floor(median(as.numeric(as.matrix(orig_seg_anal_row$Components))))
        lasso_seg_anal_row=nFactors::nScree(lasso_near_eigenvalues)
        lasso_seg_max_row=floor(median(as.numeric(as.matrix(lasso_seg_anal_row$Components))))
        lasso_seg_anal_col=nFactors::nScree(lasso_near_eigenvalues)
        lasso_seg_max_col=floor(median(as.numeric(as.matrix(lasso_seg_anal_row$Components))))
        if(verbose){message_parallel("Optimal component estimates (row):")}
        if(verbose){message_parallel(orig_seg_anal_row$Components)}
        if(verbose){message_parallel("Optimal component estimates (col):")}
        if(verbose){message_parallel(orig_seg_anal_col$Components)}
        if(verbose){message_parallel("2D LASSO Optimal component estimates (row):")}
        if(verbose){message_parallel(lasso_seg_anal_row$Components)}
        if(verbose){message_parallel("2D LASSO Optimal component estimates (col):")}
        if(verbose){message_parallel(lasso_seg_anal_col$Components)}

      }
    }
  }
  #browser()
  if(verbose){message_parallel("segmenting")}
  if(is.null(segments_orig_obj) ){
  segments_orig_obj=get2Dsegments(genomicmatrix = orig_mat_norm_small_pooled,method=seg_method,lasso=lasso_seg,smooth = smooth)}
  
  if(class(segments_orig_obj)=="integer"){segments_orig_obj=list(breakpoints_col=segments_orig_obj,breakpoints_row=segments_orig_obj,breakpoints=segments_orig_obj)}
  segments_orig_col=segments_orig_obj$breakpoints_col
  segments_orig_row=segments_orig_obj$breakpoints_row
  
  segments_orig=segments_orig_obj$breakpoints_col
  if(verbose){message_parallel("segmenting lasso")}
  if(lasso2d){
    if(is.null(segments_lasso_obj)){  

  segments_lasso_obj=get2Dsegments(genomicmatrix = lasso_cleaned_ratio_mat,method=seg_method,lasso=)
    }
  } else {
    segments_lasso_obj=segments_orig_obj
  }
  if(class(segments_lasso_obj)=="integer"){segments_lasso_obj=list(breakpoints_col=segments_lasso_obj,breakpoints_row=segments_lasso_obj,breakpoints=segments_lasso_obj)} 
  segments_lasso_col=segments_lasso_obj$breakpoints_col
  segments_lasso_row=segments_lasso_obj$breakpoints_row
  
  segments_lasso=segments_lasso_obj$breakpoints_col
  #TODO SWITCH THESE TO COL AND ROW
  segment_vector_lasso_col=foreach(i=1:(length(segments_lasso_col)+1),.combine="c",.errorhandling="pass") %do% {
    if(i==1){
      return(rep(i,segments_lasso_col[i]))
    } 
    if(i>1&i<=length(segments_lasso_col)){
      return(rep(i,diff(segments_lasso_col)[i-1]))}
    if(i>length(segments_lasso_col)){
      
      return(rep(i,ncol(lasso_cleaned_mat)-segments_lasso_col[i-1] ))}
  }
  segment_vector_lasso_row=foreach(i=1:(length(segments_lasso_row)+1),.combine="c",.errorhandling="pass") %do% {
    if(i==1){
      return(rep(i,segments_lasso_row[i]))
    } 
    if(i>1&i<=length(segments_lasso_row)){
      return(rep(i,diff(segments_lasso_row)[i-1]))}
    if(i>length(segments_lasso_row)){
      return(rep(i,nrow(lasso_cleaned_mat)-segments_lasso_row[i-1] ))}
  }
  segment_vector_orig_col=foreach(i=1:(length(segments_orig_col)+1),.combine="c",.errorhandling="pass") %do% {
    if(i==1){
      return(rep(i,segments_orig_col[i]))
    } 
    if(i>1&i<=length(segments_orig_col)){
      return(rep(i,diff(segments_orig_col)[i-1]))}
    if(i>length(segments_orig_col)){
      return(rep(i,ncol(orig_mat_norm_small_pooled)-segments_orig_col[i-1] ))}
  }
  segment_vector_orig_row=foreach(i=1:(length(segments_orig_row)+1),.combine="c",.errorhandling="pass") %do% {
    if(i==1){
      return(rep(i,segments_orig_row[i]))
    } 
    if(i>1&i<=length(segments_orig_row)){
      return(rep(i,diff(segments_orig_row)[i-1]))}
    if(i>length(segments_orig_row)){
      return(rep(i,nrow(orig_mat_norm_small_pooled)-segments_orig_row[i-1] ))}
  }
  #browser()
  orig_segment_names_col=foreach(seg_ctrl_i=1:length(unique(segment_vector_orig_col)),.combine = "c",.errorhandling="pass") %do%{
    if(seg_ctrl_i==1){(colnames(orig_mat_norm_small_pooled))[segment_vector_orig_col==seg_ctrl_i] %>% na.omit() %>% as.character() %>% CNVScope::underscored_pos_to_GRanges(zeroToOneBasedStart = F)  %>% encompass() %>% CNVScope::GRanges_to_underscored_pos()}
    if(seg_ctrl_i>1){(colnames(orig_mat_norm_small_pooled))[segment_vector_orig_col==seg_ctrl_i] %>% na.omit() %>% as.character() %>% CNVScope::underscored_pos_to_GRanges(zeroToOneBasedStart = T)  %>% encompass() %>% CNVScope::GRanges_to_underscored_pos()}
  }
  seg_ctrl_i=1
  orig_segment_names_row=foreach(seg_ctrl_i=1:length(unique(segment_vector_orig_row)),.combine = "c",.errorhandling="pass") %do%{
    if(seg_ctrl_i==1) {(rownames(orig_mat_norm_small_pooled))[segment_vector_orig_row==seg_ctrl_i] %>% na.omit() %>% as.character() %>% CNVScope::underscored_pos_to_GRanges(zeroToOneBasedStart = F)  %>% encompass() %>% CNVScope::GRanges_to_underscored_pos()}
    if(seg_ctrl_i>1) {(rownames(orig_mat_norm_small_pooled))[segment_vector_orig_row==seg_ctrl_i] %>% na.omit() %>% as.character() %>% CNVScope::underscored_pos_to_GRanges(zeroToOneBasedStart = T)  %>% encompass() %>% CNVScope::GRanges_to_underscored_pos()}
    }
  
  lasso_segment_names_col=foreach(seg_trt_i=1:length(unique(segment_vector_lasso_col)),.combine = "c",.errorhandling="pass") %do%{
    (colnames(lasso_cleaned_mat))[segment_vector_lasso_col==seg_trt_i] %>% na.omit() %>% as.character() %>% CNVScope::underscored_pos_to_GRanges(zeroToOneBasedStart = ifelse(seg_trt_i==1,F,T))  %>% encompass() %>% CNVScope::GRanges_to_underscored_pos()
    }
  lasso_segment_names_row=foreach(seg_trt_i=1:length(unique(segment_vector_lasso_row)),.combine = "c",.errorhandling="pass") %do%{
    (rownames(lasso_cleaned_mat))[segment_vector_lasso_row==seg_trt_i] %>% na.omit() %>% as.character() %>% CNVScope::underscored_pos_to_GRanges(zeroToOneBasedStart = ifelse(seg_trt_i==1,F,T))  %>% encompass() %>% CNVScope::GRanges_to_underscored_pos()}
  
  library(dtplyr)
  data.table::setDTthreads(n_cores)
  if(verbose){message_parallel("final data cleaning, step 1 (merging segments with original data)")}
  
  #creating dataframes with bins and segments for the original data, for row & col, for lasso and original:
  lasso_segment_df_small_col=tibble::tibble(loc=colnames(lasso_cleaned_mat),lasso_seg_num_col=segment_vector_lasso_col) %>% dplyr::inner_join(tibble::tibble(lasso_seg_name_col=lasso_segment_names_col,lasso_seg_num_col=1:length(lasso_segment_names_col))) 
  lasso_segment_df_small_row=tibble::tibble(loc=rownames(lasso_cleaned_mat),lasso_seg_num_row=segment_vector_lasso_row) %>% dplyr::inner_join(tibble::tibble(lasso_seg_name_row=lasso_segment_names_row,lasso_seg_num_row=1:length(lasso_segment_names_row)))
  #browser()
  orig_segment_df_small_col=tibble::tibble(loc=colnames(orig_mat_norm_small_pooled),orig_seg_num_col=segment_vector_orig_col) %>% dplyr::inner_join(tibble::tibble(orig_seg_name_col=orig_segment_names_col,orig_seg_num_col=1:length(orig_segment_names_col)))
  orig_segment_df_small_row=tibble::tibble(loc=rownames(orig_mat_norm_small_pooled),orig_seg_num_row=segment_vector_orig_row) %>% dplyr::inner_join(tibble::tibble(orig_seg_name_row=orig_segment_names_row,orig_seg_num_row=1:length(orig_segment_names_row)))
  #joining this to the original data
  (orig_mat_norm_small_pooled) %>% as.matrix() %>% reshape2::melt() %>% data.table::as.data.table() %>% data.table::setnames(c("Var1","Var2"),c("loc1","loc2")) %>% data.table::merge.data.table(orig_segment_df_small_col,by.x="loc1",by.y="loc") %>% data.table::merge.data.table(orig_segment_df_small_row,by.x="loc2",by.y="loc",suffixes = c("_loc1","_loc2") )->orig_mat_norm_small_pooled_melted  
  #joining lasso segments to the lasso data
  (lasso_cleaned_mat) %>% as.matrix() %>% reshape2::melt() %>% data.table::as.data.table() %>% data.table::setnames(c("Var1","Var2"),c("loc1","loc2")) %>% data.table::merge.data.table(lasso_segment_df_small_col,by.x="loc1",by.y="loc") %>% data.table::merge.data.table(lasso_segment_df_small_row,by.x="loc2",by.y="loc",suffixes = c("_loc1","_loc2") )->lasso_cleaned_mat_melted
  #browser()
  if(verbose){message_parallel('final data cleaning, step 2 (splitting locs)')}
  orig_mat_norm_small_pooled_melted[,c("chr1","start1","end1") := data.table::tstrsplit(loc1,"_",fixed=T)][,c("chr2","start2","end2") := data.table::tstrsplit(loc2,"_",fixed=T)][,c("chr1","end1","chr2","end2"):=NULL][,"treatment" := "ctrl"]->orig_mat_norm_small_pooled_melted_locsplit
  lasso_cleaned_mat_melted[,c("chr1","start1","end1") := data.table::tstrsplit(loc1,"_",fixed=T)][,c("chr2","start2","end2") := data.table::tstrsplit(loc2,"_",fixed=T)][,c("chr1","end1","chr2","end2"):=NULL][,"treatment" := "trt"]->lasso_cleaned_mat_melted_locsplit
  if(verbose){message_parallel('final data cleaning, step 3 (merging lasso and original data)')}
  #browser()
  small_pooled_full_merged=orig_mat_norm_small_pooled_melted_locsplit   %>% dplyr::inner_join(lasso_cleaned_mat_melted_locsplit,by=c("loc1","loc2","start1","start2"),suffix = c("_orig","_lasso"))  %>% dplyr::mutate(start1=as.numeric(start1),start2=as.numeric(start2))  %>% data.table::as.data.table()
  
  small_pooled_full_merged %>% dplyr::select(value_orig,value_lasso,dplyr::everything()) %>% tidyr::pivot_longer(cols = c(value_orig,value_lasso),names_to = "treatment",values_to = "value") %>% dplyr::mutate(treatment=gsub("value_","",treatment))  %>% dplyr::mutate(treatment=factor(treatment,levels = c("orig","lasso")))->small_pooled_full_merged_long
  if(verbose){message_parallel('final data cleaning, step 4 (adding total read counts and cells (segment pairs)) as factor variables ')}
  tibble::tibble(treatment=c("orig","lasso"),total_reads=c(sum(data.table::as.data.table(orig_mat_norm_small_pooled_melted_locsplit)$value), sum(data.table::as.data.table(lasso_cleaned_mat_melted)$value))) %>% dplyr::mutate(treatment=factor(treatment,levels = c("orig","lasso")))->total_reads_df
  #browser()  
  #get cell names for orig
  expand.grid(unique(small_pooled_full_merged$orig_seg_num_row),unique(small_pooled_full_merged$orig_seg_num_col)) %>%      dplyr::mutate(orig_cell_id= ((apply(expand.grid(unique(small_pooled_full_merged$orig_seg_name_row),unique(small_pooled_full_merged$orig_seg_name_col)), 1, function(x) paste(sort(x), collapse=" ")))))  %>% dplyr::rename(orig_seg_num_row=Var1,orig_seg_num_col=Var2)->orig_cell_df
  #get cell names for lasso
  expand.grid(unique(small_pooled_full_merged$lasso_seg_num_row),unique(small_pooled_full_merged$lasso_seg_num_col)) %>% dplyr::mutate(lasso_cell_id= ((apply(expand.grid(unique(small_pooled_full_merged$lasso_seg_name_row),unique(small_pooled_full_merged$lasso_seg_name_col)), 1, function(x) paste(sort(x), collapse=" "))))) %>% 
    
    
    #rename the segment numbers to var1 and var2
    dplyr::rename(lasso_seg_num_row=Var1,lasso_seg_num_col=Var2) %>%
    #join then by the segment numbers onto the long data.
    dplyr::inner_join(small_pooled_full_merged_long %>% data.table::as.data.table(),by=c("lasso_seg_num_row"="lasso_seg_num_row","lasso_seg_num_col"="lasso_seg_num_col")) ->lasso_cell_merged_df
  #lasso_cell_merged_df %>% dplyr::inner_join(total_reads_df,by="treatment") ->small_pooled_full_merged_with_seg_pair_ids_lasso
  
  lasso_cell_merged_df %>% dplyr::inner_join(orig_cell_df,by=c("orig_seg_num_row","orig_seg_num_col"))  %>% dplyr::rename(orig_cell_id=orig_cell_id,lasso_cell_id=lasso_cell_id)->lasso_orig_cell_merged_df
  lasso_orig_cell_merged_df %>% dplyr::inner_join(total_reads_df,by="treatment") ->small_pooled_full_merged_with_seg_pair_ids_lasso_orig
  attr(small_pooled_full_merged_with_seg_pair_ids_lasso_orig,"segments_orig_obj")<-segments_orig_obj
  attr(small_pooled_full_merged_with_seg_pair_ids_lasso_orig,"segments_lasso_obj")<-segments_lasso_obj
  return(small_pooled_full_merged_with_seg_pair_ids_lasso_orig)
}