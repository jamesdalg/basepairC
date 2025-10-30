#'  Cap and Threshold Matrix Values
#'
#' Eliminates outliers from the matrix, in both directions
#' @keywords threshold cap matrix
#' @import stats
#' @param mat matrix to be capped and thresholded
#' @param cap the maximum value to be allowed in the matrix
#' @param thresh the minimum value to be allowed in the matrix
#' @param quantile the quantile to be used as the cap value
#' @param zero_top logical, indicating whether to zero out top values
#' (if 0, then the cap value is used)
#' @return A cleaned matrix, after applying the cap and threshold
#' @examples
#' cap_thresh(matrix(1:100,10,10),cap=20,thresh=5)
#' @export

cap_thresh=function(mat,cap=20,thresh=5,quantile=0,zero_top=F){
  capped_matrix=mat
  if(quantile!=0){
    cap=quantile(as.numeric(mat[mat!=0]),quantile)
  }
  capped_matrix[abs(capped_matrix)<thresh]<-0
  if(!zero_top){
    capped_matrix[capped_matrix>0&(capped_matrix)>(cap)]<-cap
    capped_matrix[capped_matrix<0&(capped_matrix)<(-cap)]<-(-1*cap)}
  if(zero_top){
    capped_matrix[capped_matrix>0&(capped_matrix)>(cap)]<-0
    capped_matrix[capped_matrix<0&(capped_matrix)<(-cap)]<-(-1*0)
  }
  return(capped_matrix)  
}