#' Read Replicates and Process Data
#'
#' This function reads replicate files, processes the data through various steps including normalization,
#' cap and threshold adjustment, optional SCN correction, and returns either a DGEList object or a list
#' containing the processed sample matrix, conditions, and other metadata. It supports handling of Hi-C data
#' by converting genomic ranges into matrices and applying several preprocessing steps.
#' @importFrom reticulate import
#' @importFrom IRanges IRanges
#' @param replicate_files A character vector of paths to replicate files.
#' @param conditions A vector specifying the condition/group for each replicate.
#' @param start_row The starting row position for genomic range (default: 61982445).
#' @param end_row The ending row position for genomic range (default: 61988445).
#' @param start_col The starting column position for genomic range (default: 61982445).
#' @param end_col The ending column position for genomic range (default: 61988445).
#' @param bin_size The size of the bin for genomic ranges (default: 1).
#' @param chr The chromosome to be analyzed (default: "chr15").
#' @param cap The maximum value to cap the data at (default: 30).
#' @param thresh The threshold below which values are set to zero (default: 2).
#' @param quantile The quantile to use for thresholding (default: 0.9999).
#' @param return_dgel Logical; if TRUE, returns a DGEList object instead of a list (default: FALSE).
#' @param scn Logical; if TRUE, applies SCN correction to the data (default: TRUE).
#' @param make_symmetric Logical; if TRUE, makes the matrix symmetric (default: TRUE). Only applies to square matrices.
#' @param equalize_sample_coverage Logical; if TRUE, rebalances the sample coverage (default: FALSE).
#' @param n_resamples The number of resamples to use for equalizing sample coverage (default: 2).
#' @param verbose Logical; if TRUE, prints messages during processing (default: FALSE).
#' @param python Logical; if TRUE, uses Python to read the files (default: FALSE).
#' Test smaller amounts first.
#'  100 resamples will use ~600GB RAM with 112 cores and 12 samples.
#' @return If `return_dgel` is TRUE, returns a DGEList object with counts and totals.
#' Otherwise, returns a list containing the processed sample matrix, conditions, loc1, loc2,
#' underscored_positions_col, underscored_positions_row, and rep_dim.
#'
#' @examples
#' replicate_files <- c("replicate1.txt", "replicate2.txt")
#' conditions <- c("condition1", "condition2")
#' result <- read_replicates(replicate_files, conditions)
#'
#' @export
read_replicates <- function(replicate_files, conditions, start_row = 61982445, end_row = 61988445, start_col = 61982445, end_col = 61988445, bin_size = 1, chr = "chr15", cap = 30, thresh = 2, quantile = .9999, return_dgel = F, scn = T,make_symmetric=T,cores=1,reorder=T,apply_fun=pbmcapply::pbmclapply,n_resamples=2,equalize_sample_coverage=F,verbose=F,python=F) { # ,map_gc_dt=NULL
  # 

  if(reorder){
    conditions=ifelse(grepl("WT",replicate_files),"ctrl","trt")
    sample_dt=data.table::data.table(file=replicate_files,condition=conditions) %>% dplyr::arrange(condition)
    sample_dt$condition ->conditions
    sample_dt$file->replicate_files
  }
  GRanges(seqnames = chr, IRanges::IRanges(start = seq(from = start_row, to = end_row, by = bin_size), end = seq(from = start_row + bin_size - 1, to = end_row + bin_size - 1, by = bin_size))) -> row_gr
  GRanges(seqnames = chr, IRanges::IRanges(start = seq(from = start_col, to = end_col, by = bin_size), end = seq(from = start_col + bin_size - 1, to = end_col + bin_size - 1, by = bin_size))) -> col_gr
  underscored_positions_col <- CNVScope::GRanges_to_underscored_pos(col_gr)
  underscored_positions_row <- CNVScope::GRanges_to_underscored_pos(row_gr)
  if(!python){rep_dim <- data.table::fread(replicate_files[1], nThread = parallel::detectCores(), showProgress = T) %>% 
    as.matrix() %>%
    dim()} else {
      np <- reticulate::import("numpy")
      rep_dim <- np$loadtxt(replicate_files[1]) %>% 
        as.matrix() %>%
        dim()
    }
  #,fill=T
  i <- 10
  #  
  #browser()  
  apply_fun(1:length(replicate_files), function(i) {
    if(quantile==Inf|quantile=="none"){
      if(!python){
        message_parallel::message_parallel(paste0("PYTHON reading replicate ", i, " of ", length(replicate_files)))
      data.table::fread(replicate_files[i]) %>%
        as.matrix()  -> count_mat } else {
          np$loadtxt(replicate_files[i]) %>%
            as.matrix()  -> count_mat
        }
      
    } else {
      if(verbose){message_parallel(paste0("capping and thresholding replicate ", i, " of ", length(replicate_files)))}
      if(!python){
        message_parallel::message_parallel(paste0("PYTHON reading replicate ", i, " of ", length(replicate_files)))
        data.table::fread(replicate_files[i]) %>%
          as.matrix()  -> count_mat } else {
            np$loadtxt(replicate_files[i]) %>%
              as.matrix()  -> count_mat
          }
      count_mat %>%
        as.matrix() %>%
        cap_thresh(thresh = thresh, cap = cap, quantile = quantile) -> count_mat #
    }
    # data.table::fread(replicate_files[i],fill=T) ->count_dt
    # count_dt[is.na(count_dt)]<-0
    # count_dt %>%   as.matrix() %>%
    #   cap_thresh(thresh = thresh, cap = cap, quantile = quantile) -> count_mat #
    
    #count_mat[is.na(count_mat)]<-0
    
    if(verbose){basepairC::message_parallel(paste0("Read replicate ", i, " of ", length(replicate_files)))}
    if(verbose){basepairC::message_parallel(paste0("Read replicate ", replicate_files[i] ))}
    original_count_mat <- count_mat
    if (nrow(count_mat) == ncol(count_mat)) {
      if(make_symmetric){
        count_mat %>% basepairC:::make_symmetric() -> count_mat
      }
    }
    # gc & mappability correction:
    # if(!is.null(map_gc_dt)){
    if(verbose){basepairC::message_parallel("setting column and rownames on original matrix")}
    
    rownames(count_mat) <- underscored_positions_row
    colnames(count_mat) <- underscored_positions_col
    #     basepairC::message_parallel("Correcting for GC bias and mappability")
    
    # correct_mappability_gc_bias(count_mat,map_gc_dt)->count_mat
    #   }
    # SCN START
    # remove zero rows and columns
    if (scn) {
      basepairC::message_parallel("finding zero rows")
      which(Rfast::rowsums(count_mat) == 0) -> zero_rows
      which(Rfast::colsums(count_mat) == 0) -> zero_cols
      basepairC::message_parallel("setting all_row_nums")
      1:nrow(count_mat) -> all_row_nums
      1:ncol(count_mat) -> all_col_nums
      if (length(zero_rows) == 0) {
        non_zero_rows <- all_row_nums
      } else {
        all_row_nums[-zero_rows] -> non_zero_rows
      }
      if (length(zero_cols) == 0) {
        non_zero_cols <- all_col_nums
      } else {
        all_col_nums[-zero_cols] -> non_zero_cols
      }
      if(verbose){basepairC::message_parallel("subsetting count mat to remove zero rows and columns")}
      # all_row_nums[-zero_rows]->non_zero_rows
      # all_col_nums[-zero_cols]->non_zero_cols
      count_mat_without_zeros <- count_mat[non_zero_rows, non_zero_cols]
      # count mat dim:
      basepairC::message_parallel(paste0("count mat dim ", dim(count_mat)))
      basepairC::message_parallel(paste0("count mat without zeros dim ", dim(count_mat_without_zeros)))
      # run scn
      scn_count_mat <- SCN(count_mat_without_zeros + 1)
      basepairC::message_parallel(paste0("scn count mat dim ", dim(scn_count_mat)))
      # revert to original scale:
      
      scn_count_mat <- minmaxnorm(scn_count_mat) * (max(original_count_mat) - min(nonzero(original_count_mat))) + min(nonzero(original_count_mat))
      
      # add the zero rows back in, in the correct order.
      #    rbind(scn_count_mat,matrix(0,nrow = length(zero_rows),ncol=ncol(scn_count_mat)) ) -> scn_count_mat_with_zeros
      # add the zero columns back in, in the correct order.
      #    cbind(scn_count_mat_with_zeros,matrix(0,nrow = nrow(scn_count_mat_with_zeros),ncol=length(zero_cols)) ) -> scn_count_mat_with_zeros
      # reorder the rows and columns to their original order.
      # create an empty matrix of the same dimensions as the original count mat
      matrix(0, nrow = nrow(count_mat), ncol = ncol(count_mat)) -> scn_count_mat_with_zeros
      basepairC::message_parallel(paste0("scn_count_mat_with_zeros dim ", dim(scn_count_mat_with_zeros)))
      # fill in the non-zero rows and columns with the scn count mat
      scn_count_mat_with_zeros[non_zero_rows, non_zero_cols] <- scn_count_mat
      count_mat <- scn_count_mat_with_zeros
      basepairC::message_parallel(paste0("scn_count_mat_with_zeros dim ", dim(scn_count_mat_with_zeros)))
      basepairC::message_parallel(paste0("final count_mat dim ", dim(count_mat)))
      basepairC::message_parallel(paste0("length of rownames of original_count_mat ", rownames(original_count_mat)))
      basepairC::message_parallel(paste0("length of colnames of original_count_mat ", colnames(original_count_mat)))
      rownames(count_mat) <- rownames(original_count_mat)
      colnames(count_mat) <- colnames(original_count_mat)
    }
    # SCN END
    if (i == 1) {
      if(verbose){basepairC::message_parallel(paste0("after SCN dim(count_mat)", dim(count_mat)))}
      if(verbose){basepairC::message_parallel(paste0("length of  underscored_positions_col", length(underscored_positions_col)))}
      if(verbose){colnames(count_mat) <- underscored_positions_col}
      if(verbose){basepairC::message_parallel(paste0("length of  underscored_positions_row", length(underscored_positions_row)))}
      rownames(count_mat) <- underscored_positions_row
    }
    reshape2::melt(count_mat) -> count_mat_melted
    # count_mat_melted=parallel_assign_rownames(count_mat_melted,count_mat_melted$Var1,count_mat_melted$Var2)
    # return(count_mat_melted)
    options(mc.cores=cores)
    if (i != 1) {
      count_mat_melted$Var1 <- NULL
      count_mat_melted$Var2 <- NULL
      if(verbose){message_parallel(paste0("i!=1 after removing Var1 and Var2 dim(count_mat_melted)", dim(count_mat_melted)))}
      # message_parallel(paste0("stringr::str_split(basename(replicate_files[i]),pattern="-")[[1]][1]",
      #                       stringr::str_split(basename(replicate_files[i]),pattern="-")[[1]][1]
      #                        ))
      stringr::str_split(basename(replicate_files[i]), pattern = "-")[[1]][1] -> colnames(count_mat_melted)
    } else {
      if(verbose){message_parallel(paste0("i==1 dim(count_mat_melted)", dim(count_mat_melted)))}
      c("loc1", "loc2", stringr::str_split(basename(replicate_files[i]), pattern = "-")[[1]][1]) -> colnames(count_mat_melted)
    }
    return(count_mat_melted)
  }) %>% do.call(cbind, .) -> sample_mat # ,mc.cores=min(length(replicate_files),parallel::detectCores())
  #  
  if(equalize_sample_coverage){
    
    #sample_mat %>% data.table::as.data.table() %>% equalize_sample_coverage(n_samples = n_resamples,max_cores=cores) -> sample_mat
    rebalance_observed(sample_mat,filter_type = "sorted") -> sample_mat_rebalanced
    output_sample_mat=sample_mat_rebalanced
  } else {
    output_sample_mat=sample_mat
  }
  if (return_dgel) {
    DGEList(counts = as.matrix(
      data.table::as.data.table(output_sample_mat)[, !c("loc1", "loc2"), with = F]
    ), group = conditions) -> dgel
    dgel$totals <- Rfast::colsums(data.table::as.data.table(output_sample_mat)[, !c("loc1", "loc2"), with = F])
    return(dgel)
  }
  #  
  return(list(sample_mat = data.table::as.data.table(output_sample_mat)[, !c("loc1", "loc2")], conditions = conditions, loc1 = output_sample_mat$loc1, loc2 = output_sample_mat$loc2, underscored_positions_col = underscored_positions_col, underscored_positions_row = underscored_positions_row, rep_dim = rep_dim))
}