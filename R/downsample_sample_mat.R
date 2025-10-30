

#' Downsample a column of read counts.
#'
#' This function downsamples a column of read counts by retaining a specified proportion
#' of the original counts. The downsampling is done using a binomial distribution.
#'
#' @param col A numeric vector or matrix column representing the read counts to be downsampled.
#' @param prop_to_retain A numeric value or vector representing the proportion of reads to retain for each element of \code{col}. Should be between 0 and 1.
#'
#' @return A matrix with the same dimensions as \code{col} but with downsampled read counts.
#' @export
downsample_col_multinomial=function(col,min_total)
{
  #downsampled_counts <- matrix(rbinom(length(col), size = col, prob = prop_to_retain), nrow = length(col), ncol = 1)
  #downsampled_counts <- matrix(qbinom(runif(length(col)), size = col, prob = prop_to_retain), nrow = length(col), ncol = 1) #https://stackoverflow.com/questions/26241135/possible-bug-in-rbinom-for-large-numbers-of-trials
  total_counts <- sum(col)
  prob <- as.numeric(col) / total_counts
  downsampled_counts <- rmultinom(n = 1, size = min_total, prob = prob)
  colnames(downsampled_counts) <- colnames(col)
  rownames(downsampled_counts) <- rownames(col)
  return(downsampled_counts)
}

#' Downsample read counts by sample.
#'
#' This function downsamples the read counts in a sample matrix to match the sample with the lowest total read count. The downsampling is performed on a per-sample basis.
#'
#' @param sample_mat A matrix or data frame where each column represents a sample and each row represents a feature. The elements are read counts.
#' @param type A character string indicating the downsampling method. Currently, only "by_sample" is supported, which downsamples each sample to the minimum total read count across all samples.
#'
#' @return A matrix with downsampled read counts. The dimensions and names are the same as the input \code{sample_mat}.
#' @export
downsample_reads_by_sample=function(sample_mat,type="by_sample")
{
  #browser()
  total_reads=colSums(sample_mat)
  min_total=min(total_reads)
  if(!any(class(sample_mat) %in% "matrix")){sample_mat=data.frame(sample_mat)}
  #calculate the factor to downsample by
  #calculate the proportion of reads to retain-- e.g. if it's 50% greater than the minimum, then retain 66% of the reads.
  #this is to ensure that the minimum number of reads is retained
  
  if(type=="by_sample"){
    prop_to_retain=(min_total/total_reads)
    #downsampled_mat=pbapply::pblapply(1:ncol(sample_mat),function(i) {downsample_col_multinomial(sample_mat[,i],prop_to_retain[i])}) %>% do.call(cbind,.)
    downsampled_mat=pbapply::pblapply(1:ncol(sample_mat),function(i) {downsample_col_multinomial(sample_mat[,i],min_total = min_total)}) %>% do.call(cbind,.)
    downsampled_mat[is.na(downsampled_mat)]<-0
    colnames(downsampled_mat)=colnames(sample_mat)
    rownames(downsampled_mat)=rownames(sample_mat)
    
  }
  return(downsampled_mat)
}

#' Downsample counts in a DGEList object.
#'
#' This function downsamples the counts in a DGEList object to match the sample with the lowest total read count. The function works on the \code{counts} element of the DGEList.
#'
#' @param dgel A DGEList object from the \code{edgeR} package, containing the read counts to be downsampled.
#'
#' @return A DGEList object with downsampled counts.
#' @export
downsample_dgel=function(dgel){
  dgel$counts->sample_mat
  downsample_reads_by_sample(sample_mat,type="by_sample")->downsampled_mat
  dgel_downsampled=dgel
  dgel_downsampled$counts=downsampled_mat
  return(dgel_downsampled)
}
#' Downsample counts from the output of a read_replicates call.
#'
#' This function downsamples the counts in a read_replicates object to match the sample with the lowest total read count. The function works on the \code{counts} element of the read_replicates output.
#'
#' @param dgel read_replicates output list
#'
#' @return A read_replicates output list with downsampled counts.
#' @export
downsample_rr=function(rrobj){
  rrobj$sample_mat->sample_mat
  downsample_reads_by_sample(sample_mat,type="by_sample")->downsampled_mat
  rrobj_downsampled=rrobj
  rrobj_downsampled$sample_mat=downsampled_mat
  return(rrobj_downsampled)
}

downsample_col_binomial=function(col,prop_to_retain)
{
  #downsampled_counts <- matrix(rbinom(length(col), size = col, prob = prop_to_retain), nrow = length(col), ncol = 1)
  downsampled_counts <- matrix(qbinom(runif(length(col)), size = col, prob = prop_to_retain), nrow = length(col), ncol = 1) #https://stackoverflow.com/questions/26241135/possible-bug-in-rbinom-for-large-numbers-of-trials
  colnames(downsampled_counts) <- colnames(col)
  rownames(downsampled_counts) <- rownames(col)
  return(downsampled_counts)
}