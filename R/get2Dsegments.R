#' @title Get 2D Segments from a Genomic Matrix
#' @description Identifies breakpoints in a genomic matrix using segmentation algorithms, with options for smoothing and applying fused Lasso regression.
#' @param genomicmatrix A numeric matrix representing genomic data. Default is \code{NULL}.
#' @param transpose Logical; if \code{TRUE}, computes breakpoints on the transposed matrix as well. Default is \code{TRUE}.
#' @param nb_change_max Integer specifying the maximum number of breakpoints to detect. Default is \code{10}.
#' @param method Character string specifying the method used in the segmentation algorithm (e.g., "RBS"). Default is \code{"RBS"}.
#' @param smooth Logical; if \code{TRUE}, applies 2D smoothing to the genomic matrix. Default is \code{FALSE}.
#' @param smooth_window Integer specifying the window size for smoothing if \code{smooth = TRUE}. Default is \code{5}.
#' @param lasso Logical; if \code{TRUE}, applies fused Lasso regression to the genomic matrix. Default is \code{FALSE}.
#' @param lambda Numeric value for the lambda parameter in fused Lasso regression. Default is \code{NULL}.
#' @param eps Numeric value specifying the convergence threshold for fused Lasso regression. Default is \code{1e-2}.
#' @param algorithm Character string specifying the segmentation algorithm to use. Currently supports \code{"jointSeg"}. Default is \code{"jointSeg"}.
#' @param maxIter Integer specifying the maximum number of iterations for fused Lasso regression. Default is \code{20}.
#' @return A list containing the detected breakpoints:
#' \item{breakpoints_col}{Breakpoints detected in the columns of the genomic matrix.}
#' \item{breakpoints_row}{Breakpoints detected in the rows of the genomic matrix.}
#' \item{t_breakpoints_col}{(If \code{transpose = TRUE}) Breakpoints detected in the columns of the transposed genomic matrix.}
#' \item{t_breakpoints_row}{(If \code{transpose = TRUE}) Breakpoints detected in the rows of the transposed genomic matrix.}
#' @details
#' The function performs segmentation on a genomic matrix to identify breakpoints, which may indicate significant changes or patterns in the data. It allows optional smoothing and application of fused Lasso regression before segmentation.
#' \enumerate{
#'   \item If \code{smooth = TRUE}, applies 2D smoothing using \code{kernel2dsmooth} from the \code{smoothie} package on the log2-transformed genomic matrix.
#'   \item If \code{lasso = TRUE}, applies fused Lasso regression using \code{fusedLattice} from the \code{glmgen} package.
#'   \item Performs segmentation using the specified algorithm (currently only \code{"jointSeg"}) to detect breakpoints in both rows and columns.
#'   \item If \code{transpose = TRUE}, the same process is applied to the transposed genomic matrix.
#' }
#' @examples
#' \dontrun{
#' # Load necessary libraries
#' library(jointseg)
#' library(glmgen)
#' library(smoothie)
#' library(magrittr)
#'
#' # Generate a random genomic matrix
#' genomicmatrix <- matrix(rnorm(100), nrow = 10)
#'
#' # Get 2D segments without smoothing or Lasso
#' segments <- get2Dsegments(genomicmatrix)
#'
#' # Get 2D segments with smoothing
#' segments_smooth <- get2Dsegments(genomicmatrix, smooth = TRUE, smooth_window = 5)
#'
#' # Get 2D segments with Lasso regression
#' segments_lasso <- get2Dsegments(genomicmatrix, lasso = TRUE, lambda = 0.1)
#' }
#' @importFrom jointseg jointSeg
#' @importFrom glmgen fusedLattice
#' @importFrom smoothie kernel2dsmooth
#' @importFrom magrittr %>%
#' @export
get2Dsegments<-function(genomicmatrix=NULL,transpose=T,nb_change_max=10,method="RBS",smooth=F,smooth_window=5,lasso=F,lambda=NULL,eps=1e-2,algorithm="jointSeg",maxIter=20,log_before_segments=F)
{
  if(log_before_segments){
    genomicmatrix %>% log()->genomicmatrix_log;
    genomicmatrix_log[is.infinite(genomicmatrix_log)]<-0;
    genomicmatrix<-genomicmatrix_log
  }
  if(smooth){
    smoothie::kernel2dsmooth(genomicmatrix,kernel.type = "boxcar",n=smooth_window)->genomicmatrix
    
  }
  if(lasso){
    as.matrix(matrix(glmgen::fusedLattice(as.matrix(genomicmatrix), lambda = lambda, eps=eps, maxIter=maxIter)$beta,nrow=nrow(genomicmatrix))*genomicmatrix)->genomicmatrix
  }
  if(algorithm=="jointSeg"){
    #browser()
    breakpoints_col<-jointseg::jointSeg(genomicmatrix,K=nb_change_max,method=method)$bestBkp
    breakpoints_row<-jointseg::jointSeg(genomicmatrix,K=nb_change_max,method=method)$bestBkp
    if(transpose) {
      t_breakpoints_col<-jointseg::jointSeg(t(genomicmatrix),K=nb_change_max,method=method)$bestBkp
      t_breakpoints_row<-jointseg::jointSeg(t(genomicmatrix),K=nb_change_max,method=method)$bestBkp
      output_list<-list(breakpoints_col,breakpoints_row,t_breakpoints_col,t_breakpoints_row)
      names(output_list)<-c("breakpoints_col","breakpoints_row","t_breakpoints_col","t_breakpoints_row")
    }
    
  }
  return(output_list)
}
  