#' Get Segments For a Two-Condition Comparison
#'
#' This function computes the difference between two conditions and returns segmented breakpoints based on the difference.
#'
#' @param sample_mat_list A list of matrices representing sample data.
#' @param conditions A character vector indicating conditions (e.g., "ctrl" or "trt"). It must have exactly two unique values.
#' @param comp_type A character string indicating the type of comparison (options are "absdiff", "loss", "gain", or "ratio").
#'
#' @return A list containing segmented breakpoints for rows and columns, and their transposed versions.
#' @export
get_diff_segments=function(sample_mat_list,conditions=NULL,comp_type=c("absdiff","loss","gain","ratio"),return_df=F,lambda=NULL,upperQ=0.05,lowerQ=0.95){
  if(is.null(conditions)){
    tryCatch({sapply(sample_mat_list,function(x){attr(x,"sample_name")})->sample_names
      conditions=ifelse(grepl("WT",replicate_files),"ctrl","trt")},error=function(e){stop("conditions is null")})
    }
  if(length(unique(conditions))>2){stop("conditions must be length 2")}
  if(comp_type=="loss"|comp_type=="reduction"){
    Reduce("+", sample_mat_list[which("ctrl"==conditions)])-Reduce("+", sample_mat_list[which("trt"==conditions)])->diff_mat #possibly use DescTools::Windsorise
    #browser()
    cap=quantile(neg_vals(nonzero(as.numeric(diff_mat))),upperQ)
    thresh=quantile(neg_vals(nonzero(as.numeric(diff_mat))),lowerQ)
    
  }

  if(comp_type=="gain"|comp_type=="increase"){
    Reduce("+", sample_mat_list[which("trt"==conditions)])-Reduce("+", sample_mat_list[which("ctrl"==conditions)])->diff_mat
    cap=quantile(nonzero(as.numeric(diff_mat)),upperQ)
    thresh=quantile(nonzero(as.numeric(diff_mat)),lowerQ)
    
  }
  if(comp_type=="absdiff"){
    abs(Reduce("+", sample_mat_list[which("trt"==conditions)])-Reduce("+", sample_mat_list[which("ctrl"==conditions)]))->diff_mat
    cap=quantile(nonzero(as.numeric(diff_mat)),upperQ)
    thresh=quantile(nonzero(as.numeric(diff_mat)),lowerQ)
    
  }
  if(comp_type=="ratio"){
    create_ratio_mat(Reduce("+", sample_mat_list[which("ctrl"==conditions)]),Reduce("+", sample_mat_list[which("trt"==conditions)]),cap=Inf,thresh=0)->diff_mat
    cap=quantile(nonzero(as.numeric(diff_mat)),upperQ)
    thresh=quantile(nonzero(as.numeric(diff_mat)),lowerQ)
    
  }
    
    preprocess_mat(diff_mat,verbose=T,dsfactor=1,lambda=NULL,cap=cap,thresh=thresh)->diff_mat_df
  seg_obj=attr(diff_mat_df,"segments_orig_obj")
  if(return_df){attr(seg_obj,"diff_mat_df")=diff_mat_df}
  return(seg_obj)
}
neg_vals=function(x){
  x[x<0]
}
pos_vals=function(x){
  x[x>0]
}
signed_log_transform=function(x){
  sign(x)*log(abs(x)+1)
}