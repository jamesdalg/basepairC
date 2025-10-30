#' Create Long-Form Control and Treatment Data Frame
#'
#' This function reads, preprocesses, and segments control and treatment matrices
#' to create a long-form data frame. It performs various data processing steps
#' and segmentation to prepare the data for further analysis.
#'
#' @param ctl_filename The filename of the control matrix.
#' @param trt_filename The filename of the treatment matrix.
#' @param bin_size The bin size for genomic binning (default is 3).
#' @param thresh The threshold value for data processing (default is 5).
#' @param cap The cap value for data processing (default is 20).
#' @param lambda The regularization parameter for lasso filtering (default is 5).
#' @param eps The convergence threshold for lasso filtering (default is 1e-10).
#' @param maxIter The maximum number of iterations for lasso filtering (default is 1e4).
#' @param quantile The quantile value for data processing (default is 0.99).
#' @param zero_top Logical, indicating whether to zero out top values (default is TRUE).
#' @param chr The chromosome name (default is "chr15").
#' @param start_row The start position of the genomic region (default is 61982445).
#' @param end_row The end position of the genomic region (default is 61988445).
#' @param start_col The start position of the genomic region (default is 61982445).
#' @param end_col The end position of the genomic region (default is 61988445).
#' @param dsfactor The downsampling factor (default is 1).
#' @param mat_start The start index for matrix subsetting (default is 1).
#' @param mat_end The end index for matrix subsetting (default is 2000).
#' @param plot Logical, indicating whether to generate plots (default is FALSE).
#' @param ctrl_seg_max The maximum number of segments for control (default is 9).
#' @param trt_seg_max The maximum number of segments for treatment (default is 7).
#' @param lasso_ratio_seg_max The maximum number of segments for lasso ratio (default is 7).
#' @param ctrl_seg_method The segmentation method for control (default is "B").
#' @param trt_seg_method The segmentation method for treatment (default is "B").
#' @param lasso_ratio_seg_method The segmentation method for lasso ratio (default is "G").
#' @param ctrl_seg_model The segmentation model for control (default is "Dplus").
#' @param trt_seg_model The segmentation model for treatment (default is "Dplus").
#' @param lasso_ratio_seg_model The segmentation model for lasso ratio (default is "Dplus").
#' @param debug Logical, indicating whether to enter debug mode (default is FALSE).
#' @param make_symmetric Logical, indicating whether to make the matrix symmetric (default is TRUE).
#' @param seg_algorithm The segmentation algorithm to use (default is "HiCseg"). jointSeg also works
#' @param ctrl_lasso Perform 2D LASSO on ctrl?
#' @param trt_lasso Perform 2D LASSO on trt? 
#' @param total_read_correction Perform total read correction?
#' @param parallel_message Logical, indicating whether to display parallel message (default is TRUE).
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom magrittr %>%
#' @importFrom data.table := tstrsplit
#' @importFrom smoothie kernel2dsmooth
#' @return A processed and segmented data frame.
#'
#' @examples
#' ctl_file <- "control_matrix.csv"
#' trt_file <- "treatment_matrix.csv"
#' df <- create_long_ctrl_trt_df(ctl_file, trt_file)
#'
#' @export
create_long_ctrl_trt_df_sym=function(ctl_filename,trt_filename,bin_size=3,thresh=5,cap=20,lambda=5,eps=1e-10,maxIter=20,quantile=0,zero_top=T,
                                 chr="chr15",start_row=61982445,end_row=61988445,start_col=61982445,end_col=61988445,dsfactor=1,mat_start_row=1,mat_end_row=2000,mat_start_col=1,
                                 mat_end_col=2000,plot=F,ctrl_seg_max=9,trt_seg_max=7,lasso_ratio_seg_max=7,ctrl_seg_method="G",trt_seg_method="G",lasso_ratio_seg_method="G",
                                 ctrl_seg_model="Dplus",trt_seg_model="Dplus",lasso_ratio_seg_model="Dplus",debug=F,n_cores=1,make_symmetric=T,ds_fn="sum",seg_algorithm="HiCseg",ctrl_lasso=T,trt_lasso=T,total_read_correction=F,parallel_message=T){
  #browser()
  
  options(nThread=n_cores)
  if(debug){browser()}
  print('reading data')
  if(parallel_message){message_parallel("Reading data")}
  if(any("matrix" %in% class(ctl_filename) )){
    print("reading matrix data, enforcing symmetry")
    if(make_symmetric){
      gene_ctrl_mat_norm=ctl_filename %>% make_symmetric() %>% as("sparseMatrix") %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col] %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)
      gene_trt_mat_norm=trt_filename %>% make_symmetric() %>% as("sparseMatrix") %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col] %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)} else{
        gene_ctrl_mat_norm=ctl_filename  %>% as("sparseMatrix") %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col] %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)
        gene_trt_mat_norm=trt_filename  %>% as("sparseMatrix") %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col]  %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)
      }
  }
  else{
    if(make_symmetric){
      #browser()
      print("reading matrix data")
      gene_ctrl_mat_norm=(data.table::fread(ctl_filename,nThread=n_cores,showProgress = T) %>% as.matrix() %>% make_symmetric() %>% as("sparseMatrix")) %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col] %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)
      gene_trt_mat_norm=(data.table::fread(trt_filename,nThread=n_cores,showProgress = T) %>% as.matrix() %>% make_symmetric() %>% as("sparseMatrix")) %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col] %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)} else {
        gene_ctrl_mat_norm=(data.table::fread(ctl_filename,nThread=n_cores,showProgress = T) %>% as.matrix()  %>% as("sparseMatrix")) %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col] %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)
        gene_trt_mat_norm=(data.table::fread(trt_filename,nThread=n_cores,showProgress = T) %>% as.matrix()  %>% as("sparseMatrix")) %>% .[mat_start_row:mat_end_row,mat_start_col:mat_end_col]  %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)
      }
  }
  
  thresh=thresh*dsfactor
  # if(ctrl_lasso){
  #   #browser()
  #   print("performing 2D lasso on control matrix")
  #   if(parallel_message){message_parallel("Performing 2D lasso on control matrix")}
  #   as.matrix(matrix(glmgen::fusedLattice(as.matrix(gene_ctrl_mat_norm), lambda = lambda, eps=eps, maxIter=maxIter)$beta,nrow=nrow(gene_ctrl_mat_norm))*gene_ctrl_mat_norm)->gene_ctrl_mat_norm  
  # }
  # if(trt_lasso){
  #   #browser()
  #   print("performing 2D lasso on treatment matrix")
  #   if(parallel_message){message_parallel("Performing 2D lasso on treatment matrix")}
  #   as.matrix(matrix(glmgen::fusedLattice(as.matrix(gene_trt_mat_norm), lambda = lambda, eps=eps, maxIter=maxIter)$beta,nrow=nrow(gene_trt_mat_norm))*gene_trt_mat_norm)->gene_trt_mat_norm
  # }
  #  bin_size=3
  GRanges(seqnames=chr, IRanges(start=seq(from=start_row,to=end_row,by=bin_size), end=seq(from=start_row+bin_size-1,to=end_row+bin_size-1,by=bin_size)))->row_gr
  GRanges(seqnames=chr, IRanges(start=seq(from=start_col,to=end_col,by=bin_size), end=seq(from=start_col+bin_size-1,to=end_col+bin_size-1,by=bin_size)))->col_gr
  underscored_positions_col=CNVScope::GRanges_to_underscored_pos(col_gr)
  underscored_positions_row=CNVScope::GRanges_to_underscored_pos(row_gr)
  #browser()
  (row_gr %>% CNVScope::GRanges_to_underscored_pos())[1:dim(gene_ctrl_mat_norm)[1]]->rownames(gene_ctrl_mat_norm)
  (col_gr %>% CNVScope::GRanges_to_underscored_pos())[1:dim(gene_ctrl_mat_norm)[2]]->colnames(gene_ctrl_mat_norm)
  (row_gr %>% CNVScope::GRanges_to_underscored_pos())[1:dim(gene_trt_mat_norm)[1]]->rownames(gene_trt_mat_norm)
  (col_gr %>% CNVScope::GRanges_to_underscored_pos())[1:dim(gene_trt_mat_norm)[2]]->colnames(gene_trt_mat_norm)
  #col_gr %>% CNVScope::GRanges_to_underscored_pos()->colnames(gene_trt_mat_norm)->rownames(gene_trt_mat_norm)
  #dsfactor=1
  n_cores=parallel::detectCores()/2
  print('downsampling')
  if(parallel_message){message_parallel("Downsampling")}
  #browser()
  CNVScope::downsample_genomic_matrix(as.matrix(gene_ctrl_mat_norm),downsamplefactor = dsfactor,singlechromosome = T)->gene_ctrl_mat_norm_small
  #gene_ctrl_mat_norm_small_pooled <- sum_pooling_parallel(gene_ctrl_mat_norm, factor=dsfactor, ncores=n_cores) %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)
  if(ds_fn=="sum"){gene_ctrl_mat_norm_small_pooled <- block_sum_matrix(gene_ctrl_mat_norm, dsfactor, dsfactor) %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)} 
  if(ds_fn=="mean"){gene_ctrl_mat_norm_small_pooled <- gene_ctrl_mat_norm_small %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)}
  colnames(gene_ctrl_mat_norm_small)->colnames(gene_ctrl_mat_norm_small_pooled)
  rownames(gene_ctrl_mat_norm_small)->rownames(gene_ctrl_mat_norm_small_pooled)
  CNVScope::downsample_genomic_matrix(as.matrix(gene_trt_mat_norm),downsamplefactor = dsfactor,singlechromosome = T)->gene_trt_mat_norm_small
  #gene_trt_mat_norm_small_pooled <- sum_pooling_parallel(gene_trt_mat_norm, factor=dsfactor,ncores=n_cores) %>% cap_thresh(thresh=thresh,quantile=quantile)
  if(ds_fn=="sum"){gene_trt_mat_norm_small_pooled <- block_sum_matrix(gene_trt_mat_norm, dsfactor, dsfactor) %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)}
  if(ds_fn=="mean"){gene_trt_mat_norm_small_pooled <- gene_trt_mat_norm_small %>% cap_thresh(thresh=thresh, quantile=quantile, cap=cap)}
  #browser()
  colnames(gene_trt_mat_norm_small)->colnames(gene_trt_mat_norm_small_pooled)
  rownames(gene_trt_mat_norm_small)->rownames(gene_trt_mat_norm_small_pooled)
  if(ctrl_lasso){
    #browser()
    print("performing 2D lasso on control matrix")
    if(parallel_message){message_parallel("Performing 2D lasso on control matrix")}
    gene_ctrl_mat_norm_small_pooled->gene_ctrl_mat_norm_small_pooled_original
    as.matrix(matrix(glmgen::fusedLattice(as.matrix(gene_ctrl_mat_norm_small_pooled), lambda = lambda, eps=eps, maxIter=maxIter)$beta,nrow=nrow(gene_ctrl_mat_norm_small_pooled))*gene_ctrl_mat_norm_small_pooled)->gene_ctrl_mat_norm_small_pooled
    gene_ctrl_mat_norm_small_pooled->gene_ctrl_mat_norm_small_pooled_lasso
  } else{gene_ctrl_mat_norm_small_pooled_original=gene_ctrl_mat_norm_small_pooled}
  if(trt_lasso){
    #browser()
    print("performing 2D lasso on treatment matrix")
    if(parallel_message){message_parallel("Performing 2D lasso on treatment matrix")}
    gene_trt_mat_norm_small_pooled->gene_trt_mat_norm_small_pooled_original
    as.matrix(matrix(glmgen::fusedLattice(as.matrix(gene_trt_mat_norm_small_pooled), lambda = lambda, eps=eps, maxIter=maxIter)$beta,nrow=nrow(gene_trt_mat_norm_small_pooled))*gene_trt_mat_norm_small_pooled)->gene_trt_mat_norm_small_pooled
    gene_trt_mat_norm_small_pooled->gene_trt_mat_norm_small_pooled_lasso
  } else{gene_trt_mat_norm_small_pooled_original=gene_trt_mat_norm_small_pooled}
  if(total_read_correction){
    #browser()
    print('total read correction')
    total_read_ratio_list=total_read_correction_list(list(gene_ctrl_mat_norm_small_pooled,gene_trt_mat_norm_small_pooled))
    gene_ctrl_mat_norm_small_pooled=total_read_ratio_list[[1]]
    gene_trt_mat_norm_small_pooled=total_read_ratio_list[[2]]
  }
  print('segmentation')
  if(parallel_message){message_parallel("Segmentation")}
  if(debug){browser()}
  if(is.null(ctrl_seg_max)){
    ctrl_near_eigenvalues=eigen(Matrix::nearPD(cov(gene_ctrl_mat_norm_small_pooled))$mat)$values
    #browser()
    ctrl_seg_anal=nFactors::nScree(ctrl_near_eigenvalues)
    ctrl_seg_max=floor(median(as.numeric(as.matrix(ctrl_seg_anal$Components))))
    print("Optimal component estimates (ctrl):")
    if(parallel_message){message_parallel("Optimal component estimates (ctrl):")}
    print(ctrl_seg_anal$Components)
    if(parallel_message){message_parallel(ctrl_seg_anal$Components)}
  }
  if(is.null(trt_seg_max)){
    trt_near_eigenvalues=eigen(Matrix::nearPD(cov(gene_trt_mat_norm_small_pooled))$mat)$values
    trt_seg_anal=nFactors::nScree(trt_near_eigenvalues)
    trt_seg_max=floor(median(as.numeric(as.matrix(trt_seg_anal$Components))))
    print("Optimal component estimates (trt):")
    print(trt_seg_anal$Components)
    if(parallel_message){message_parallel("Optimal component estimates (trt):")}
    if(parallel_message){message_parallel(trt_seg_anal$Components)}
  }
  
  segments_ctrl_obj=CNVScope::getAsymmetricBlockIndices(genomicmatrix = gene_ctrl_mat_norm_small_pooled,distrib = ctrl_seg_method,model=ctrl_seg_model,nb_change_max = ctrl_seg_max,algorithm=seg_algorithm)
  if(debug){browser()}
  if(class(segments_ctrl_obj)=="integer"){segments_ctrl_obj=list(breakpoints_col=segments_ctrl_obj,breakpoints_row=segments_ctrl_obj,breakpoints=segments_ctrl_obj)}
  segments_ctrl_col=segments_ctrl_obj$breakpoints_col
  segments_ctrl_row=segments_ctrl_obj$breakpoints_row
  if(length(segments_ctrl_col)==0){segments_ctrl_col=c(1,ncol(gene_ctrl_mat_norm_small_pooled))} else {segments_ctrl_col=segments_ctrl_obj$breakpoints_col}
  if(length(segments_ctrl_row)==0){segments_ctrl_row=c(1,nrow(gene_ctrl_mat_norm_small_pooled))} else {segments_ctrl_row=segments_ctrl_obj$breakpoints_row}
  segments_ctrl=segments_ctrl_col #segments_ctrl_obj$breakpoints_col
  segments_trt_obj=CNVScope::getAsymmetricBlockIndices(genomicmatrix = gene_trt_mat_norm_small_pooled,distrib = trt_seg_method,model=trt_seg_model,nb_change_max = trt_seg_max,algorithm = seg_algorithm)
  #browser()  
  if(class(segments_trt_obj)=="integer"){segments_trt_obj=list(breakpoints_col=segments_trt_obj,breakpoints_row=segments_trt_obj)}
  segments_trt_col=segments_trt_obj$breakpoints_col
  segments_trt_row=segments_trt_obj$breakpoints_row
  if(length(segments_trt_col)==0){segments_trt_col=c(1,ncol(gene_trt_mat_norm_small_pooled))} else {segments_trt_col=segments_trt_obj$breakpoints_col}
  if(length(segments_trt_row)==0){segments_trt_row=c(1,nrow(gene_trt_mat_norm_small_pooled))} else {segments_trt_row=segments_trt_obj$breakpoints_row}
  
  segments_trt=segments_trt_row #segments_trt_obj$breakpoints_row
  #browser()
  # segment_vector_trt=foreach::foreach(i=1:length(segments_trt),.combine="c") %do% {
  #   if(i==1){
  #     rep(i,segments_trt[i])
  #   } else{ rep(i,diff(segments_trt)[i-1])}
  # }
  # segment_vector_ctrl=foreach::foreach(i=1:length(segments_ctrl),.combine="c") %do% {
  #   if(i==1){
  #     rep(i,segments_ctrl[i])
  #   } else{ rep(i,diff(segments_ctrl)[i-1])}
  #   
  # }
  segment_vector_trt=foreach(i=1:(length(segments_trt)+1),.combine="c") %do% {
    #browser()
    if(i==1){
      #browser()
      return(rep(i,segments_trt[i]))
    } 
    if(i>1&i<=length(segments_trt)){
      #browser()
      return(rep(i,diff(segments_trt)[i-1]))}
    if(i>length(segments_trt)){
      
      return(rep(i,ncol(gene_trt_mat_norm_small_pooled)-segments_trt[i-1] ))}
  }
  segment_vector_ctrl=foreach(i=1:(length(segments_ctrl)+1),.combine="c") %do% {
    #browser()
    if(i==1){
      #browser()
      #print(rep(i,segments_ctrl[i]))
      #print(segments_ctrl[i])
      return(rep(i,segments_ctrl[i]))
    } 
    if(i>1&i<=length(segments_ctrl)){
      #browser()
      #print(rep(i,diff(segments_ctrl)[i-1]))
      #print(diff(segments_ctrl)[i-1])
      return(rep(i,diff(segments_ctrl)[i-1]))}
    if(i>length(segments_ctrl)){
      #print(rep(i,ncol(gene_ctrl_mat_norm_small_pooled)-segments_ctrl[i-1] ))
      #print(ncol(gene_ctrl_mat_norm_small_pooled)-segments_ctrl[i-1])
      return(rep(i,ncol(gene_ctrl_mat_norm_small_pooled)-segments_ctrl[i-1] ))}
  }
  ctrl_segment_names=foreach(seg_ctrl_i=1:length(unique(segments_ctrl)),.combine = "c") %do%{
    (colnames(gene_ctrl_mat_norm_small_pooled))[segment_vector_ctrl==seg_ctrl_i] %>% na.omit() %>% as.character() %>% CNVScope::underscored_pos_to_GRanges()  %>% encompass() %>% CNVScope::GRanges_to_underscored_pos()}
  
  trt_segment_names=foreach(seg_trt_i=1:length(unique(segments_trt)),.combine = "c") %do%{
    (colnames(gene_trt_mat_norm_small_pooled))[segment_vector_trt==seg_trt_i] %>% na.omit() %>% as.character() %>% CNVScope::underscored_pos_to_GRanges()  %>% encompass() %>% CNVScope::GRanges_to_underscored_pos()}
  print('creating ratio matrix and lasso filtering')
  if(parallel_message){message_parallel('creating ratio matrix and lasso filtering')}
  #this creates the ratio matrix and does the lasso filtering with lambda=5 and epsilon 1e-10. Lambda is the most potent parameter of the two to change.
  trt_ratio_mat=create_ratio_mat(gene_ctrl_mat_norm_small_pooled,gene_trt_mat_norm_small_pooled,quantile=quantile,thresh=thresh,cap=cap,zero_top = T)
  #TODO: add documentation using formulae with direct multiplication symbol and lambda/epsilon.
  #browser()
  if(is.null(lambda)){
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
  }
  as.matrix(matrix(glmgen::fusedLattice(as.matrix(trt_ratio_mat), lambda = lambda, eps=eps, maxIter=maxIter)$beta,nrow=nrow(trt_ratio_mat))*trt_ratio_mat)->lasso_cleaned_ratio_mat
  
  if(debug){browser()}
  #if(is.null(lasso_ratio_seg_max)){lasso_ratio_seg_max=nFactors::nScree(lasso_cleaned_ratio_mat)}
  colnames(gene_trt_mat_norm_small)->colnames(lasso_cleaned_ratio_mat)
  rownames(gene_trt_mat_norm_small)->rownames(lasso_cleaned_ratio_mat)
  #this actually does the segmentation. Feel free to adjust nb_change_max as seems prudent. If the segments don't cover homogeneous areas, increase it. If the segments have tiny ones driven by outliers, cap it more thoroughly (e.g. set quantile=0.99 instead of .99) and decrease nb_change_max once outliers are adequately controlled.
  if(is.null(lasso_ratio_seg_max)){
    #browser()
    eigen(Matrix::nearPD(lasso_cleaned_ratio_mat)$mat)$values->lr_eigs
    lr_lasso_seg_anal=nFactors::nScree(lr_eigs)
    if(plot){plot(1:length(lr_eigs),lr_eigs,xlab="Eigenvalue Number",ylab="Eigenvalue",main="Scree Plot of Eigenvalues for 2D graphical lasso ratio(ctrl/trt) \n Choose the value at the curve elbow to determine the maximum number of segmentation change points.",type="bs")}
    lasso_ratio_seg_max=floor(median(as.numeric(as.matrix(lr_lasso_seg_anal$Components))))
    print("Optimal component estimates (2D graphical lasso ratio(ctrl/trt)):")
    print(lr_lasso_seg_anal$Components)
  } 
  #browser()
  smooth_window=5
  lasso_cleaned_ratio_mat %>% log2()->lasso_cleaned_ratio_mat_log2;
  lasso_cleaned_ratio_mat_log2[is.infinite(lasso_cleaned_ratio_mat_log2)]<-0;
  smoothie::kernel2dsmooth(lasso_cleaned_ratio_mat_log2,kernel.type = "boxcar",n=smooth_window)->lasso_cleaned_ratio_mat_log2_smooth
  segments_trt_ratio_lasso=CNVScope::getAsymmetricBlockIndices(genomicmatrix = lasso_cleaned_ratio_mat_log2_smooth,distrib = "G",model="D",nb_change_max = lasso_ratio_seg_max,algorithm=seg_algorithm)
  #1:nrow(lasso_cleaned_ratio_mat)->seg_numbers
  #ifelse(seg_numbers %in% segments_trt_ratio_lasso$breakpoints_row,1,0)->seg_numbers
  #ComplexHeatmap::Heatmap(lasso_cleaned_ratio_mat,show_column_dend = F,show_column_names = F,show_heatmap_legend = F,show_row_names = T,cluster_rows = F,cluster_columns = F,row_labels=1:nrow(lasso_cleaned_ratio_mat))
  rm(lasso_cleaned_ratio_mat_log2);rm(lasso_cleaned_ratio_mat_log2_smooth);
  if(class(segments_trt_ratio_lasso)=="integer"){
    segments_trt_ratio_lasso=list(breakpoints_row=segments_trt_ratio_lasso,
                                  breakpoints_col=segments_trt_ratio_lasso)
  }
  segments_trt_ratio_lasso_col=segments_trt_ratio_lasso$breakpoints_col
  segments_trt_ratio_lasso_row=segments_trt_ratio_lasso$breakpoints_row
  segments_trt_ratio_lasso=segments_trt_ratio_lasso$breakpoints_col
  #this gets the indexes of the segments in the ratio lasso trt comparison.
  #browser()
  segment_vector_trt_ratio_lasso=foreach(i=1:(length(segments_trt_ratio_lasso)+1),.combine="c") %do% {
    #browser()
    if(i==1){
      #browser()
      return(rep(i,segments_trt_ratio_lasso[i]))
    } 
    if(i>1&i<=length(segments_trt_ratio_lasso)){
      #browser()
      return(rep(i,diff(segments_trt_ratio_lasso)[i-1]))}
    if(i>length(segments_trt_ratio_lasso)){
      
      return(rep(i,ncol(lasso_cleaned_ratio_mat)-segments_trt_ratio_lasso[i-1] ))}
  }
  segment_vector_trt_ratio_lasso_row=foreach(i=1:(length(segments_trt_ratio_lasso_row)+1),.combine="c") %do% {
    #browser()
    if(i==1){
      #browser()
      return(rep(i,segments_trt_ratio_lasso_row[i]))
    } 
    if(i>1&i<=length(segments_trt_ratio_lasso_row)){
      #browser()
      return(rep(i,diff(segments_trt_ratio_lasso_row)[i-1]))}
    if(i>length(segments_trt_ratio_lasso_row)){
      
      return(rep(i,nrow(lasso_cleaned_ratio_mat)-segments_trt_ratio_lasso_row[i-1] ))}
  }
  segment_vector_trt_ratio_lasso_col=foreach(i=1:(length(segments_trt_ratio_lasso_col)+1),.combine="c") %do% {
    #browser()
    if(i==1){
      #browser()
      return(rep(i,segments_trt_ratio_lasso_col[i]))
    } 
    if(i>1&i<=length(segments_trt_ratio_lasso_col)){
      #browser()
      return(rep(i,diff(segments_trt_ratio_lasso_col)[i-1]))}
    if(i>length(segments_trt_ratio_lasso_col)){
      
      return(rep(i,ncol(lasso_cleaned_ratio_mat)-segments_trt_ratio_lasso_col[i-1] ))}
  }
  #browser()
  #this gets the unique segment underscored positions as a character vector.
  trt_ratio_lasso_segment_names=foreach(seg_trt_rl_i=1:length(unique(segment_vector_trt_ratio_lasso)),.combine = "c") %do%{
    
    #browser()
    print(seg_trt_rl_i)
    (colnames(gene_trt_mat_norm_small_pooled))[segment_vector_trt_ratio_lasso==seg_trt_rl_i] %>% na.omit() %>% as.character() %>% CNVScope::underscored_pos_to_GRanges()  %>% encompass() %>% CNVScope::GRanges_to_underscored_pos()}
  if(plot==T){
    
    #plot lasso_cleaned_ratio_mat with ComplexHeatmap
    print("plotting lasso_cleaned_ratio_mat")
    #browser()
    top_anno=ComplexHeatmap::HeatmapAnnotation(df=data.frame(segment_lasso_ratio=letters[segment_vector_trt_ratio_lasso_col]),which="column")
    side_anno=ComplexHeatmap::rowAnnotation(df=data.frame(segment_lasso_ratio=letters[segment_vector_trt_ratio_lasso_row]))
    #lasso_hm_list=ComplexHeatmap::Heatmap(CNVScope::signedRescale(as.matrix(lasso_cleaned_ratio_mat),max_cap = 25),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="ctrl/trt",top_annotation = top_anno,bottom_annotation = top_anno,right_annotation = side_anno,left_annotation = side_anno) 
    #as.matrix(1/lasso_cleaned_ratio_mat)->lasso_cleaned_ratio_mat_recip
    #lasso_cleaned_ratio_mat_recip[is.infinite(lasso_cleaned_ratio_mat_recip)]<-1
    #lasso_hm_list=ComplexHeatmap::Heatmap(CNVScope::signedRescale(as.matrix(lasso_cleaned_ratio_mat),max_cap = 25),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="trt/ctrl",top_annotation = top_anno,bottom_annotation = top_anno,right_annotation = side_anno,left_annotation = side_anno)
    lasso_hm_list=ComplexHeatmap::Heatmap(CNVScope::signedRescale(as.matrix(lasso_cleaned_ratio_mat),max_cap = 25),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="trt/ctrl",top_annotation = top_anno,bottom_annotation = top_anno,right_annotation = side_anno,left_annotation = side_anno)
    #+ ComplexHeatmap::HeatmapAnnotation(df=data.frame(segment_trt=letters[segment_vector_trt_ratio_lasso]),which="row")
    #+ ComplexHeatmap::HeatmapAnnotation(df=data.frame(segment_trt=letters[segment_vector_trt_ratio_lasso]),which="column")
    png("lasso_cleaned_ratio_mat.png",width=2000,height=2000,res=300)
    ComplexHeatmap::draw(lasso_hm_list)
    dev.off()
    #png("lasso_cleaned_ratio_mat_recip.png",width=2000,height=2000,res=300)
    #ComplexHeatmap::draw(lasso_hm_list_recip)
    #dev.off()
    print("plotting trt_ratio_mat")
    #ComplexHeatmap::draw(lasso_hm_list)
    print("plotting original data")
    #browser()
    # ComplexHeatmap::Heatmap(CNVScope::signedRescale(as.matrix(gene_ctrl_mat_norm),max_cap = 25),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="ctrl (original)") %>% ComplexHeatmap::draw()
    # ComplexHeatmap::Heatmap(CNVScope::signedRescale(log2(as.matrix(gene_ctrl_mat_norm)+1),max_cap = 25),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="log2(ctrl+1) (original)") %>% ComplexHeatmap::draw()
    ComplexHeatmap::Heatmap(as.matrix(gene_ctrl_mat_norm),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="ctrl (original)") %>% ComplexHeatmap::draw()
    ComplexHeatmap::Heatmap(log2(as.matrix(gene_ctrl_mat_norm)+1),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="log2(ctrl+1) (original)") %>% ComplexHeatmap::draw()
    print("orig summary:")
    print(summary(as.numeric(gene_ctrl_mat_norm_small_pooled_original)))
    if(ctrl_lasso){
      print("orig lasso summary")
      print(summary(as.numeric(gene_ctrl_mat_norm_small_pooled_lasso))) 
    }
    ComplexHeatmap::Heatmap(as.matrix(gene_trt_mat_norm),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="trt (original)") %>% ComplexHeatmap::draw()
    ComplexHeatmap::Heatmap(log2(as.matrix(gene_trt_mat_norm)+1),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="log2(trt+1) (original)") %>% ComplexHeatmap::draw()
    # ComplexHeatmap::Heatmap(CNVScope::signedRescale(as.matrix(gene_trt_mat_norm),max_cap = 25),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="trt (original)") %>% ComplexHeatmap::draw()
    # ComplexHeatmap::Heatmap(CNVScope::signedRescale(log2(as.matrix(gene_trt_mat_norm)+1),max_cap = 25),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="log2(trt+1) (original)") %>% ComplexHeatmap::draw()
    print("trt summary:")
    print(summary(as.numeric(gene_trt_mat_norm_small_pooled_original)))
    if(trt_lasso){
      print("trt lasso summary")
      print(summary(as.numeric(gene_trt_mat_norm_small_pooled_lasso)))}
    # ComplexHeatmap::Heatmap(CNVScope::signedRescale(as.matrix(gene_ctrl_mat_norm_small_pooled),max_cap = 25),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="ctrl (lasso)",top_annotation = top_anno,bottom_annotation = top_anno,right_annotation = side_anno,left_annotation = side_anno)  %>% ComplexHeatmap::draw()
    # ComplexHeatmap::Heatmap(CNVScope::signedRescale(log2(as.matrix(gene_ctrl_mat_norm_small_pooled)+1),max_cap = 25),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="log2(ctrl) (lasso)",top_annotation = top_anno,bottom_annotation = top_anno,right_annotation = side_anno,left_annotation = side_anno)  %>% ComplexHeatmap::draw()
    ComplexHeatmap::Heatmap(as.matrix(gene_ctrl_mat_norm_small_pooled),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="ctrl (lasso)",top_annotation = top_anno,bottom_annotation = top_anno,right_annotation = side_anno,left_annotation = side_anno)  %>% ComplexHeatmap::draw()
    ComplexHeatmap::Heatmap(log2(as.matrix(gene_ctrl_mat_norm_small_pooled)+1),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="log2(ctrl) (lasso)",top_annotation = top_anno,bottom_annotation = top_anno,right_annotation = side_anno,left_annotation = side_anno)  %>% ComplexHeatmap::draw()
    
    # ComplexHeatmap::Heatmap(CNVScope::signedRescale(as.matrix(gene_trt_mat_norm_small_pooled),max_cap = 25),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="trt (lasso)",top_annotation = top_anno,bottom_annotation = top_anno,right_annotation = side_anno,left_annotation = side_anno) %>% ComplexHeatmap::draw()
    # ComplexHeatmap::Heatmap(CNVScope::signedRescale(log2(as.matrix(gene_trt_mat_norm_small_pooled)+1),max_cap = 25),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="log2(trt+1) (lasso)",top_annotation = top_anno,bottom_annotation = top_anno,right_annotation = side_anno,left_annotation = side_anno) %>% ComplexHeatmap::draw()
    ComplexHeatmap::Heatmap(as.matrix(gene_trt_mat_norm_small_pooled),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="trt (lasso)",top_annotation = top_anno,bottom_annotation = top_anno,right_annotation = side_anno,left_annotation = side_anno) %>% ComplexHeatmap::draw()
    ComplexHeatmap::Heatmap(log2(as.matrix(gene_trt_mat_norm_small_pooled)+1),cluster_rows = F,cluster_columns = F,show_row_names = F,show_column_names = F,column_title="log2(trt+1) (lasso)",top_annotation = top_anno,bottom_annotation = top_anno,right_annotation = side_anno,left_annotation = side_anno) %>% ComplexHeatmap::draw()
    
    #ComplexHeatmap::draw(lasso_hm_list_recip)
  }
  library(dtplyr)
  data.table::setDTthreads(parallel::detectCores())
  print("final data cleaning, step 1")
  if(parallel_message){message_parallel("final data cleaning, step 1")}
  #creating dataframes with bins and segments for the original data: trt, trt_ratio_lasso, and ctrl
  #browser()
  trt_segment_lasso_ratio_df_small=tibble::tibble(loc=colnames(lasso_cleaned_ratio_mat),trt_ratio_lasso_seg_num=segment_vector_trt_ratio_lasso) %>%   dplyr::inner_join(tibble::tibble(trt_ratio_lasso_seg_name=trt_ratio_lasso_segment_names,trt_ratio_lasso_seg_num=1:length(trt_ratio_lasso_segment_names)))
  trt_segment_df_small=tibble::tibble(loc=colnames(gene_trt_mat_norm_small_pooled),trt_seg_num=segment_vector_trt) %>% dplyr::inner_join(tibble::tibble(trt_seg_name=trt_segment_names,trt_seg_num=1:length(trt_segment_names))) 
  ctrl_segment_df_small=tibble::tibble(loc=colnames(gene_ctrl_mat_norm_small_pooled),ctrl_seg_num=segment_vector_ctrl) %>% dplyr::inner_join(tibble::tibble(ctrl_seg_name=ctrl_segment_names,ctrl_seg_num=1:length(ctrl_segment_names)))
  #joining this to the original data
  #ctrl
  #(gene_ctrl_mat_norm_small_pooled) %>% reshape2::melt() %>% dplyr::rename(loc1=Var1,loc2=Var2) %>% dtplyr::lazy_dt() %>% dplyr::inner_join(ctrl_segment_df_small,by=c("loc1"="loc")) %>% dplyr::inner_join(ctrl_segment_df_small,by=c("loc2"="loc"),suffix =  c("_loc1","_loc2") )->gene_ctrl_mat_norm_small_pooled_melted
  #data.table::as.data.table(gene_ctrl_mat_norm_small_pooled) %>% data.table::melt.data.table(id.vars="loc1",variable.name = 'loc2')
  if(debug){browser()}
  #browser()
  (gene_ctrl_mat_norm_small_pooled) %>% as.matrix() %>% reshape2::melt() %>% data.table::as.data.table() %>% data.table::setnames(c("Var1","Var2"),c("loc1","loc2")) %>% data.table::merge.data.table(ctrl_segment_df_small,by.x="loc1",by.y="loc") %>% data.table::merge.data.table(ctrl_segment_df_small,by.x="loc2",by.y="loc",suffixes = c("_loc1","_loc2") )->gene_ctrl_mat_norm_small_pooled_melted
  #trt
  #(gene_trt_mat_norm_small_pooled) %>% reshape2::melt() %>% dplyr::rename(loc1=Var1,loc2=Var2) %>% dtplyr::lazy_dt() %>% dplyr::inner_join(trt_segment_df_small,by=c("loc1"="loc")) %>% dplyr::inner_join(trt_segment_df_small,by=c("loc2"="loc"),suffix =  c("_loc1","_loc2") )->gene_trt_mat_norm_small_pooled_melted
  (gene_trt_mat_norm_small_pooled) %>% as.matrix() %>% reshape2::melt() %>% data.table::as.data.table() %>% data.table::setnames(c("Var1","Var2"),c("loc1","loc2")) %>% data.table::merge.data.table(trt_segment_df_small,by.x="loc1",by.y="loc") %>% data.table::merge.data.table(trt_segment_df_small,by.x="loc2",by.y="loc",suffixes = c("_loc1","_loc2") )->gene_trt_mat_norm_small_pooled_melted
  #trt Ratio
  #(lasso_cleaned_ratio_mat) %>% reshape2::melt() %>% dplyr::rename(loc1=Var1,loc2=Var2) %>% dtplyr::lazy_dt() %>% dplyr::inner_join(trt_segment_lasso_ratio_df_small,by=c("loc1"="loc")) %>% dplyr::inner_join(trt_segment_lasso_ratio_df_small,by=c("loc2"="loc"),suffix =  c("_loc1","_loc2") )->gene_trt_lr_mat_norm_small_pooled_melted
  (lasso_cleaned_ratio_mat) %>% reshape2::melt() %>% data.table::as.data.table() %>% data.table::setnames(c("Var1","Var2"),c("loc1","loc2")) %>% data.table::merge.data.table(trt_segment_lasso_ratio_df_small,by.x="loc1",by.y="loc") %>% data.table::merge.data.table(trt_segment_lasso_ratio_df_small,by.x="loc2",by.y="loc",suffixes = c("_loc1","_loc2") )->gene_trt_lr_mat_norm_small_pooled_melted
  print("final data cleaning, step 2")
  if(parallel_message){message_parallel("final data cleaning, step 2")}
  #separating out chr1,start1,end1,chr2,start2,end2 for all three matrices from chr_1_100 chr_100_200, etc.
  
  #gene_ctrl_mat_norm_small_pooled_melted  %>% tidyr::separate(loc1,into=c("chr1","start1","end1"),sep = "_",remove = F) %>% dplyr::select(-end1) %>% dplyr::select(-chr1) %>% tidyr::separate(loc2,into=c("chr2","start2","end2"),sep = "_",remove = F) %>% dplyr::select(-end2) %>% dplyr::select(-chr2) -> gene_ctrl_mat_norm_small_pooled_melted_merged
  #gene_trt_mat_norm_small_pooled_melted  %>% tidyr::separate(loc1,into=c("chr1","start1","end1"),sep = "_",remove = F) %>% dplyr::select(-end1) %>% dplyr::select(-chr1) %>% tidyr::separate(loc2,into=c("chr2","start2","end2"),sep = "_",remove = F) %>% dplyr::select(-end2) %>% dplyr::select(-chr2) -> gene_trt_mat_norm_small_pooled_melted_merged
  #gene_trt_lr_mat_norm_small_pooled_melted  %>% tidyr::separate(loc1,into=c("chr1","start1","end1"),remove = F,sep = "_") %>% dplyr::select(-chr1) %>% dplyr::select(-end1) %>% tidyr::separate(loc2,into=c("chr2","start2","end2"),remove = F,sep = "_") %>% dplyr::select(-chr2) %>% dplyr::select(-end2)  %>% dplyr::filter(value!=0) %>% na.omit() -> gene_trt_lr_mat_norm_small_pooled_melted_merged
  #rewrite in data.table
  #browser()
  gene_ctrl_mat_norm_small_pooled_melted[,c("chr1","start1","end1") := tstrsplit(loc1,"_",fixed=T)][,c("chr2","start2","end2") := tstrsplit(loc2,"_",fixed=T)][,c("chr1","end1","chr2","end2"):=NULL][,"treatment" := "ctrl"]->gene_ctrl_mat_norm_small_pooled_melted_merged
  gene_trt_mat_norm_small_pooled_melted[,c("chr1","start1","end1") := tstrsplit(loc1,"_",fixed=T)][,c("chr2","start2","end2") := tstrsplit(loc2,"_",fixed=T)][,c("chr1","end1","chr2","end2"):=NULL][,"treatment" := "trt"]->gene_trt_mat_norm_small_pooled_melted_merged
  gene_trt_lr_mat_norm_small_pooled_melted[,c("chr1","start1","end1") := tstrsplit(loc1,"_",fixed=T)][,c("chr2","start2","end2") := tstrsplit(loc2,"_",fixed=T)][,c("chr1","end1","chr2","end2"):=NULL][,"treatment" := "trtlr"]->gene_trt_lr_mat_norm_small_pooled_melted_merged
  
  #merging the three matrices together
  print("final data cleaning, step 3")
  if(parallel_message){message_parallel("final data cleaning, step 3")}
  small_pooled_full_merged=gene_ctrl_mat_norm_small_pooled_melted_merged   %>% dplyr::inner_join(gene_trt_mat_norm_small_pooled_melted_merged,by=c("loc1","loc2","start1","start2"),suffix = c("_ctrl","_trt"))  %>% dplyr::mutate(start1=as.numeric(start1),start2=as.numeric(start2)) %>%  dplyr::inner_join(gene_trt_lr_mat_norm_small_pooled_melted_merged %>% dplyr::mutate(start1=as.numeric(start1),start2=as.numeric(start2)),by=c("loc1","loc2"),suffix = c("","_trtlr")) %>% dplyr::rename(value_trtlr=value) %>% data.table::as.data.table()
  #rewrite in data.table
  #gene_ctrl_mat_norm_small_pooled_melted_merged %>% data.table::as.data.table() %>% data.table::merge.data.table(gene_trt_mat_norm_small_pooled_melted_merged,by=c("loc1","loc2","start1","start2"),suffixes = c("_ctrl","_trt"))  %>% data.table::merge.data.table(gene_trt_lr_mat_norm_small_pooled_melted_merged %>% data.table::setnames(c("start1","start2"),c("start1_trtlr","start2_trtlr")),by=c("loc1","loc2"),suffixes = c("","_trtlr")) %>% data.table::setnames(c("start1_trtlr","start2_trtlr"),c("start1","start2")) ->small_pooled_full_merged_dt
  #making a long form.
  small_pooled_full_merged %>% dplyr::select(-treatment) %>% dplyr::select(value_ctrl,value_trt,value_trtlr,dplyr::everything()) %>% tidyr::pivot_longer(cols = c(value_ctrl,value_trt,value_trtlr),names_to = "treatment",values_to = "value") %>% dplyr::mutate(treatment=gsub("value_","",treatment))  %>% dplyr::mutate(treatment=factor(treatment,levels = c("ctrl","trt","trtlr")))->small_pooled_full_merged_long
  #adding cell id and then getting the total read counts for each treatment.
  options(copy=T)
  print("final data cleaning, step 4")
  if(parallel_message){message_parallel("final data cleaning, step 4")}
  if(debug){browser()}
  tibble::tibble(treatment=c("ctrl","trt","trtlr"),total_reads=c(sum(data.table::as.data.table(gene_ctrl_mat_norm_small_pooled_melted_merged)$value), sum(data.table::as.data.table(gene_trt_mat_norm_small_pooled_melted_merged)$value),sum(data.table::as.data.table(gene_trt_lr_mat_norm_small_pooled_melted_merged)$value))) %>% dplyr::mutate(treatment=factor(treatment,levels = c("ctrl","trt","trtlr")))->total_reads_df
  #tibble::tibble(treatment=c("ctrl","trt","trtlr"),total_reads=c(sum(gene_ctrl_mat_norm_small_pooled), sum(gene_trt_mat_norm_small_pooled),sum(lasso_cleaned_ratio_mat))) %>% dplyr::mutate(treatment=factor(treatment,levels = c("ctrl","trt","trtlr")))->total_reads_df
  
  expand.grid(unique(small_pooled_full_merged$trt_ratio_lasso_seg_num_loc1),unique(small_pooled_full_merged$trt_ratio_lasso_seg_num_loc2)) %>% 
    dplyr::mutate(trtlr_cell_id= ((apply(expand.grid(unique(small_pooled_full_merged$trt_ratio_lasso_seg_name_loc1),unique(small_pooled_full_merged$trt_ratio_lasso_seg_name_loc2)), 1, function(x) paste(sort(x), collapse=" "))))) %>% 
    dplyr::rename(trtlr_seg_num_loc1=Var1,trtlr_seg_num_loc2=Var2) %>% dplyr::inner_join(small_pooled_full_merged_long %>% data.table::as.data.table(),by=c("trtlr_seg_num_loc1"="trt_ratio_lasso_seg_num_loc1","trtlr_seg_num_loc2"="trt_ratio_lasso_seg_num_loc2")) ->trt_cell_merged_df
  trt_cell_merged_df %>% dplyr::inner_join(total_reads_df,by="treatment") ->small_pooled_full_merged_with_seg_pair_ids_trtlr
  return(small_pooled_full_merged_with_seg_pair_ids_trtlr)
}