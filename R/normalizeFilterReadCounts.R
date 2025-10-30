#' Normalize and Filter Read Counts
#'
#' This function normalizes read counts from 2D basepairC data, optionally filters the data,
#' and returns various outputs based on the specified output type. It supports returning
#' raw DGEList objects, normalized counts, counts per million (CPM), and log-transformed CPM.
#' It can also split the normalized counts back into matrices for multiple conditions.
#'
#' @param sample_mat A matrix of raw read counts where rows represent genes and columns represent samples.
#' @param conditions A factor or vector that specifies the condition/group for each column in `sample_mat`.
#' @param loc1 A vector indicating the first location identifier for each row in `sample_mat`.
#' @param loc2 A vector indicating the second location identifier for each row in `sample_mat`.
#' @param nrow_output The desired number of rows in the output matrices.
#' @param ncol_output The desired number of columns in the output matrices.
#' @param filter Logical; if TRUE, performs filtering based on median counts per million.
#' @param underscored_positions_col Column names for the output matrices, typically indicating positions with underscores.
#' @param underscored_positions_row Row names for the output matrices, typically indicating positions with underscores.
#' @param output_type A character string specifying the type of output: "dgel" for DGEList object, "corrected_counts" for normalized counts,
#' "cpm" for counts per million, or "logcpm" for log-transformed counts per million. Defaults to "corrected_counts".
#' @param norm_factor_type A character string specifying the type of normalization factor to use: "TMMwsp" for TMM with singleton pairing,
#'  "TMM" for TMM normalization, or "DESeq2" for DESeq2 normalization, or "TRC" for total read correction. Defaults to "TMMwsp", which works well for zero inflated data.
#'  "NC" for normalizing by the number of nonzero counts.
#' @param filter_type A character string specifying the type of filtering to use: "ave_ab" for average log counts per million, "mean_loc_sum" for mean of the sum of counts per location,
#' or "diff" for the difference between the sum of counts in one condition versus all other conditions. Defaults to "ave_ab".
#' @param estimate_dispersion Logical; if TRUE, estimates dispersion using edgeR's `estimateDisp` function.
#' @param norm_offsets Logical; if TRUE, calculates normalization offsets using csaw's `normOffsets` function.
#' @param skip_mats Logical; if TRUE, skips the creation of matrices and returns the raw data table with indices.
#' @importFrom csaw normFactors normOffsets 
#' @importFrom edgeR estimateDisp
#'
#' @return Depending on the `output_type` and `filter` parameters, this function can return a DGEList object,
#' a list of matrices with sum and variance for each condition, or a data table with normalized (filtered) counts.
#' For "dgel" output type, returns a DGEList object. For "corrected_counts", "cpm", or "logcpm" with `filter=FALSE`,
#' returns a list containing matrices of sum and variance for each condition, along with the conditions vector.
#' If `filter=TRUE`, returns a list containing matrices of sum and variance for filtered data, conditions vector,
#' and a data table linking locations to indices.
#'
#' @examples
#' # Assuming `sample_mat`, `conditions`, `loc1`, `loc2` are defined:
#' result <- normalizeFilterReadCounts(sample_mat, conditions, loc1, loc2, 100, 5)
#' # To get normalized counts without filtering:
#' normalized_counts <- normalizeFilterReadCounts(sample_mat, conditions, loc1, loc2, 100, 5, output_type = "corrected_counts")
#'
#' @export
normalizeFilterReadCounts <- function(sample_mat, conditions, loc1, loc2, nrow_output, ncol_output, filter = F, underscored_positions_col = underscored_positions_col, underscored_positions_row, output_type = "corrected_counts",norm_factor_type="TMMwsp",filter_type="ave_ab",estimate_dispersion=F,norm_offsets=F,skip_mats=F) {
  # convert each matrix to a numeric vector
  # Combine matrices and create group factor


  # Create DGEList object
  print("creating DGElist Object")
  dgel <- DGEList(counts = sample_mat, group = conditions)
  print("adding totals for each sample")
  dgel$totals <- Rfast::colsums(as.matrix(data.table::as.data.table(sample_mat)))
  print("calculating normalization factors for each sample")
  #csaw::normFactors(dgel) -> dgel
  if(norm_factor_type=="TMMwsp"){
    edgeR::calcNormFactors(dgel,method="TMMwsp") -> dgel
  }
  if(norm_factor_type=="TMM"){
    csaw::normFactors(dgel) -> dgel
  }
  if(norm_factor_type=="DESeq2"){
    dgel$samples$norm.factors=DESeq2::estimateSizeFactorsForMatrix(dgel$counts)
    #DESeq2::estimateSizeFactors(dgel) -> dgel
  }
  if(norm_factor_type=="TRC"){
      sums=Rfast::colsums(dgel$counts)
      dgel$samples$norm.factors=1/(sums/min(sums))
      #dgel$samples$norm.factors=
    #DESeq2::estimateSizeFactors(dgel) -> dgel
  }
  if(norm_factor_type=="NC"){
    #scaling by the number of nonzero counts.
    nz_mat=dgel$counts
    nz_mat[nz_mat!=0]<-1
    sums=Rfast::colsums(dgel$counts)
    dgel$samples$norm.factors=sums/(sums/min(sums))
  }
  
  #edgeR::calcNormFactors(dgel,method="TMMwsp") -> dgel
  print("outputting normalized counts")
  if (output_type == "dgel" && filter == F) {
    return(dgel)
  }
  if (output_type == "corrected_counts") {
    if(norm_factor_type=="DESeq2"|norm_factor_type=="TRC"|norm_factor_type=="NC"){
      #browser()
     #sweep(dgel$counts, 2, dgel$samples$norm.factors, FUN="/")->normalized_counts
      Rfast::eachrow(dgel$counts,dgel$samples$norm.factors,"/")->normalized_counts
      colnames(normalized_counts) <- colnames(dgel$counts)
    } 
    if(norm_factor_type %in% c("TMMwsp","TMM")){
    normalized_counts <- convert_to_normalized_counts(dgel)
    }
    #normalized_counts <- sweep(dgel$counts, 2, dgel$samples$norm.factors, FUN="*")
  }
  if (tolower(output_type) == "cpm") {
    normalized_counts <- cpm(dgel, log = FALSE, normalized.lib.sizes = T)
  }
  if (tolower(output_type) == "logcpm") {
    normalized_counts <- cpm(dgel, log = FALSE, normalized.lib.sizes = T)
  }


  # Reference:https://support.bioconductor.org/p/132240/

  # glmQLFit(dgel_filtered_disp)->dgel_filtered_fit

  if (filter == F) {
    # Split normalized counts back into a series of matrices (so it will work for multiple treatments at the same time).
    # browser()
    n_conditions <- length(unique(conditions))
    # options(mc.cores=min(parallel::detectCores(),n_conditions))
    print("combining counts for each condition into a matrix")
    pbapply::pblapply(1:n_conditions, function(i) {
      matrix(Rfast::rowsums(as.matrix(normalized_counts[, which(conditions == unique(conditions)[i])])), nr = nrow_output, nc = ncol_output) -> output_mat
      colnames(output_mat) <- underscored_positions_col
      rownames(output_mat) <- underscored_positions_row
      return(output_mat)
    }) -> sum_mat_list
    print("calculating variance for each condition")
    pbapply::pblapply(1:n_conditions, function(i) {
      matrix(Rfast::rowVars(as.matrix(normalized_counts[, which(conditions == unique(conditions)[i])])), nr = nrow_output, nc = ncol_output) -> output_mat
      colnames(output_mat) <- underscored_positions_col
      rownames(output_mat) <- underscored_positions_row
      return(output_mat)
    }) -> var_mat_list

    # control_normalized <- normalized_counts[, 1:ncol(control)]
    # treatment_normalized <- normalized_counts[, (ncol(control) + 1):ncol(all_counts)]
    # matrix(Rfast::rowsums(control_normalized), nrow = nrow_output, ncol = ncol(control)) -> control_normalized
    # matrix(treatment_normalized, nrow = nrow(treatment), ncol = ncol(treatment)) -> treatment_normalized

    # Return as a list of matrices
    print("Complete!")
    return(list(sum_mat_list = sum_mat_list, var_mat_list = var_mat_list, conditions = conditions))
  } else {
    print("calculating average log counts per million")
    #browser()
    #set zero to na in counts
    #count_mat=dgel$counts
    #count_mat[count_mat==0]<-NA
    
    ave.ab <- aveLogCPM(dgel, normalized.lib.sizes = F)
    # browser()
    #old way
    if(filter_type=="ave_ab"){
    ave.ab > median(ave.ab) -> keep}
    #new way
    if(filter_type=="mean_loc_sum"){
     Rfast::rowsums(dgel$counts)->sums 
      sums > mean(sums)   -> keep
    }
    #diff based
    if(filter_type=="diff"){
      #this needs to work for n unique conditions.
      #split the counts by condition
      lapply(1:length(unique(conditions)),function(i){
        dgel$counts[,which(conditions==unique(conditions)[i])]->condition_counts
        dgel$counts[,which(conditions!=unique(conditions)[i])]->non_condition_counts
        Rfast::rowsums(condition_counts)-Rfast::rowsums(non_condition_counts)->diffs
        as.numeric(diffs > 0) -> keep_i
        return(keep_i)
      }) %>% do.call(cbind,.) %>% Rfast::rowsums() ->keep_vec
      keep_vec>0->keep
      # dgel$counts[,which(conditions==unique(conditions)[i])]->condition_counts
      # dgel$counts[,which(conditions!=unique(conditions)[i])]->non_condition_counts
      # Rfast::rowsums(condition_counts)-Rfast::rowsums(non_condition_counts)->diffs
      # diffs > 0 -> keep
    }
    
    #Matrix::rowMeans(count_mat,na.rm=T)->ave.nz.count
    
    print("removing those squares less than the median CPM")
    dgel_filtered <- DGEList(counts = as.matrix(
      data.table::as.data.table(sample_mat)[keep, ]
    ), group = conditions)
    print("calculating norm factors")
    if(norm_factor_type=="TMMwsp"){
      edgeR::calcNormFactors(dgel_filtered,method="TMMwsp") -> dgel_filtered
    }
    if(norm_factor_type=="TMM"){
      csaw::normFactors(dgel_filtered) -> dgel_filtered
    }
    if(norm_factor_type=="DESeq2"){
      dgel_filtered$samples$norm.factors=DESeq2::estimateSizeFactorsForMatrix(dgel_filtered$counts)
      #DESeq2::estimateSizeFactors(dgel_filtered) -> dgel_filtered
    }
    if(norm_factor_type=="TRC"){
      #dgel_filtered$samples$norm.factors=DESeq2::estimateSizeFactorsForMatrix(dgel_filtered$counts)
      sums=Rfast::colsums(dgel_filtered$counts)
      dgel_filtered$samples$norm.factors=1/(sums/min(sums))
      #DESeq2::estimateSizeFactors(dgel_filtered) -> dgel_filtered
    }
    if(norm_factor_type=="NC"){
      #scaling by the number of nonzero counts.
      nz_mat=dgel_filtered$counts
      nz_mat[nz_mat!=0]<-1
      sums=Rfast::colsums(dgel_filtered$counts)
      dgel_filtered$samples$norm.factors=sums/(sums/min(sums))
    }
    #csaw::normFactors(dgel_filtered) -> dgel_filtered
    #edgeR::calcNormFactors(dgel_filtered,method="TMMwsp") -> dgel_filtered
    if(norm_offsets){
    print("calculating norm offsets")
    csaw::normOffsets(dgel_filtered, se.out = T) -> dgel_filtered}
    if(estimate_dispersion){
    print("estimating dispersion")
    edgeR::estimateDisp(dgel_filtered, design = model.matrix(~conditions), trend.method = "movingave") -> dgel_filtered}
    if (output_type == "dgel") {
      return(dgel_filtered)
    }
     #browser()
    # cpm(dgel_filtered,log=FALSE,normalized.lib.sizes = T)->normalized_filtered_sample_matrix
    if (output_type == "corrected_counts") {
      print("outputting normalized counts")
      if(norm_factor_type=="DESeq2"|norm_factor_type=="TRC"|norm_factor_type=="NC"){
        #sweep(dgel_filtered$counts, 2, dgel_filtered$samples$sizeFactors, FUN="/")->normalized_filtered_sample_matrix
        #Rfast::eachrow(dgel_filtered$counts,dgel_filtered$samples$norm.factors,"/")->normalized_counts
        Rfast::eachrow(dgel_filtered$counts,dgel_filtered$samples$norm.factors,"/") %>% data.table::as.data.table()->normalized_filtered_sample_matrix
        colnames(normalized_filtered_sample_matrix) <- colnames(dgel_filtered$counts)
        #browser()
      } else{
      normalized_filtered_sample_matrix <- convert_to_normalized_counts(dgel_filtered) %>% data.table::as.data.table()
      }
    }
    if (tolower(output_type) == "cpm") {
      print("outputting normalized CPM")
      normalized_filtered_sample_matrix <- cpm(dgel_filtered, log = FALSE, normalized.lib.sizes = T) %>% data.table::as.data.table()
    }
    if (tolower(output_type) == "logcpm") {
      print("outputting normalized log CPM")
      normalized_filtered_sample_matrix <- cpm(dgel_filtered, log = FALSE, normalized.lib.sizes = T) %>% data.table::as.data.table()
    }
    # normalized_filtered_sample_matrix<-data.table::as.data.table(normalized_filtered_sample_matrix)
    #browser()
    normalized_filtered_sample_matrix$loc1 <- as.character(loc1)[keep]
    normalized_filtered_sample_matrix$loc2 <- as.character(loc2)[keep]
    # create empty sample matrix and fill it with the counts at the locations indicated by keep.
    # data.table::as.data.table(normalized_counts)->empty_sample_mat #copy over the shape of the original matrix.
    empty_sample_mat <- data.table::as.data.table(matrix(0, nrow = nrow(normalized_counts), ncol = ncol(normalized_counts)))
    colnames(empty_sample_mat) <- colnames(normalized_counts)
    # empty_sample_mat[TRUE]<-0 #zero everything out.
    empty_sample_mat$loc1 <- as.character(loc1)
    empty_sample_mat$loc2 <- as.character(loc2)
    empty_sample_mat[which(keep), ] <- normalized_filtered_sample_matrix # add back in only the filtered bits.
    #browser()
    rr_obj=list(sample_mat=empty_sample_mat[,!c("loc1","loc2"),with=F],conditions=conditions,loc1=empty_sample_mat$loc1,loc2=empty_sample_mat$loc2,underscored_positions_col=underscored_positions_col,underscored_positions_row=underscored_positions_row,rep_dim=c(nrow_output,ncol_output))
    # browser()
    if(skip_mats){
      loc_index_tbl <- data.table::data.table(loc1 = loc1[keep], loc2 = loc2[keep], index = which(keep))

      return(list(sum_mat_list = NULL, var_mat_list = NULL, offset_sum_mat_list=NULL, offset_mean_mat_list=NULL,
                  mean_nonzero_mat_list = NULL, med_nonzero_mat_list = NULL, 
                  conditions = conditions, loc_index_tbl = loc_index_tbl,initial_sample_mat=sample_mat,final_sample_mat=empty_sample_mat,dgel_filtered=dgel_filtered,rr_obj=rr_obj))
    }

    n_conditions <- length(unique(conditions))
    # options(mc.cores=min(parallel::detectCores(),n_conditions))
    pbmcapply::pbmclapply(1:n_conditions, function(i) {
      matrix(Rfast::rowsums(as.matrix(
        empty_sample_mat[, colnames(empty_sample_mat)[which(conditions == unique(conditions)[i])], with = FALSE]
      )), nr = nrow_output, nc = ncol_output) -> output_mat
      colnames(output_mat) <- underscored_positions_col
      rownames(output_mat) <- underscored_positions_row
      attributes(output_mat)$condition=unique(conditions)[i]
      attributes(output_mat)$sample_names=colnames(empty_sample_mat)[which(conditions == unique(conditions)[i])]
      return(output_mat)
    }, mc.cores = min(parallel::detectCores(), n_conditions)) -> sum_mat_list_filtered
    sapply(sum_mat_list_filtered,function(mat){attributes(mat)$condition})->names(sum_mat_list_filtered)
    #names(sum_mat_list_filtered)=conditions
    #create a new matrix the same size, set zeros to NA. then take the nonzero mean of each row
    empty_sample_mat_na_nonzero=empty_sample_mat
    empty_sample_mat_na_nonzero[empty_sample_mat_na_nonzero==0]<-NA
    pbapply::pblapply(1:n_conditions, function(i) {
      matrix(Matrix::rowMeans(na.rm = T,x = as.matrix(
        empty_sample_mat_na_nonzero[, colnames(empty_sample_mat)[which(conditions == unique(conditions)[i])], with = FALSE]
      )), nr = nrow_output, nc = ncol_output) -> output_mat
      output_mat[is.na(output_mat)]<-0
      colnames(output_mat) <- underscored_positions_col
      rownames(output_mat) <- underscored_positions_row
      attributes(output_mat)$condition=unique(conditions)[i]
      attributes(output_mat)$sample_names=colnames(empty_sample_mat)[which(conditions == unique(conditions)[i])]
      return(output_mat)
    } ) -> mean_nonzero_mat_list_filtered#floor(min(parallel::detectCores(), n_conditions)/n_conditions) #, mc.cores = min(parallel::detectCores(), n_conditions)
    #names(mean_nonzero_mat_list_filtered)=conditions
    sapply(mean_nonzero_mat_list_filtered,function(mat){attributes(mat)$condition})->names(mean_nonzero_mat_list_filtered)
    #do the same for the median.
    empty_sample_mat_na_nonzero=empty_sample_mat
    empty_sample_mat_na_nonzero[empty_sample_mat_na_nonzero==0]<-NA
    pbapply::pblapply(1:n_conditions, function(i) {
      matrix(Rfast::rowMedians(na.rm = T,x = as.matrix(
        empty_sample_mat_na_nonzero[, colnames(empty_sample_mat)[which(conditions == unique(conditions)[i])], with = FALSE]
      )), nr = nrow_output, nc = ncol_output) -> output_mat
      output_mat[is.na(output_mat)]<-0
      colnames(output_mat) <- underscored_positions_col
      rownames(output_mat) <- underscored_positions_row
      attributes(output_mat)$condition=unique(conditions)[i]
      attributes(output_mat)$sample_names=colnames(empty_sample_mat)[which(conditions == unique(conditions)[i])]
      return(output_mat)
    } ) -> med_nonzero_mat_list_filtered #, mc.cores = min(parallel::detectCores(), n_conditions)
    #names(med_nonzero_mat_list_filtered)=conditions
    sapply(med_nonzero_mat_list_filtered,function(mat){attributes(mat)$condition})->names(med_nonzero_mat_list_filtered)
        pbmcapply::pbmclapply(1:n_conditions, function(i) {
      matrix(Rfast::rowVars(as.matrix(
        empty_sample_mat[, colnames(empty_sample_mat)[which(conditions == unique(conditions)[i])], with = FALSE]
      )), nr = nrow_output, nc = ncol_output) -> output_mat
      colnames(output_mat) <- underscored_positions_col
      rownames(output_mat) <- underscored_positions_row
      attributes(output_mat)$condition=unique(conditions)[i]
      attributes(output_mat)$sample_names=colnames(empty_sample_mat)[which(conditions == unique(conditions)[i])]
      return(output_mat)
    }) -> var_mat_list_filtered
    
    #names(var_mat_list_filtered)=conditions
    sapply(var_mat_list_filtered,function(mat){attributes(mat)$condition})->names(var_mat_list_filtered)
#browser()
  if(!is.null(dgel_filtered$offset)){  
    offset_dt=dgel_filtered$offset %>% data.table::as.data.table()
    pbmcapply::pbmclapply(1:n_conditions, function(i) {
      matrix(Rfast::rowsums(as.matrix(
        offset_dt[, colnames(empty_sample_mat)[which(conditions == unique(conditions)[i])], with = FALSE]
      )), nr = nrow_output, nc = ncol_output) -> output_mat
      colnames(output_mat) <- underscored_positions_col
      rownames(output_mat) <- underscored_positions_row
      attributes(output_mat)$condition=unique(conditions)[i]
      attributes(output_mat)$sample_names=colnames(empty_sample_mat)[which(conditions == unique(conditions)[i])]
      return(output_mat)
    }, mc.cores = min(parallel::detectCores(), n_conditions)) -> offset_sum_mat_list_filtered
    sapply(offset_sum_mat_list_filtered,function(mat){attributes(mat)$condition})->names(offset_sum_mat_list_filtered)
    pbmcapply::pbmclapply(1:n_conditions, function(i) {
      matrix(Rfast::rowmeans(as.matrix(
        offset_dt[, colnames(empty_sample_mat)[which(conditions == unique(conditions)[i])], with = FALSE]
      )), nr = nrow_output, nc = ncol_output) -> output_mat
      colnames(output_mat) <- underscored_positions_col
      rownames(output_mat) <- underscored_positions_row
      attributes(output_mat)$condition=unique(conditions)[i]
      attributes(output_mat)$sample_names=colnames(empty_sample_mat)[which(conditions == unique(conditions)[i])]
      return(output_mat)
    }, mc.cores = min(parallel::detectCores(), n_conditions)) -> offset_mean_mat_list_filtered
    sapply(offset_mean_mat_list_filtered,function(mat){attributes(mat)$condition})->names(offset_mean_mat_list_filtered)
  } else { offset_sum_mat_list_filtered=NULL;offset_mean_mat_list_filtered=NULL  }

  
    loc_index_tbl <- data.table::data.table(loc1 = loc1[keep], loc2 = loc2[keep], index = which(keep))
    
    # setDT(loc_index_tbl)
    # loc_index_tbl[, =: index=1:nrow(.)]
    return(list(sum_mat_list = sum_mat_list_filtered, var_mat_list = var_mat_list_filtered, offset_sum_mat_list=offset_sum_mat_list_filtered, offset_mean_mat_list=offset_mean_mat_list_filtered,
                 mean_nonzero_mat_list = mean_nonzero_mat_list_filtered, med_nonzero_mat_list = med_nonzero_mat_list_filtered, 
                conditions = conditions, loc_index_tbl = loc_index_tbl,initial_sample_mat=sample_mat,final_sample_mat=empty_sample_mat,dgel_filtered=dgel_filtered))
    # if("DGEList" %in% return_vals){
    # return(dgel_filtered)}
    # }
  }
}
