## ----setup, include=FALSE,warning=F-------------------------------------------
knitr::opts_chunk$set(echo = TRUE,cache=T,fig.width=12, fig.height=8)
library(magrittr)
library(basepairC)

## ----read_replicates,eval=F---------------------------------------------------
# list.files(system.file(package="basepairC","data/example_matrices"),pattern="*.txt",recursive=T,full.names = T)->replicate_files
# conditions=ifelse(grepl("WT",replicate_files),"ctrl","trt")
# sample_dt=data.table::data.table(file=replicate_files,condition=conditions) %>% dplyr::arrange(condition)
# sample_dt$condition ->conditions
# suppressMessages(basepairC::read_replicates(replicate_files,conditions,thresh=0,cap=Inf,scn = F,apply_fun=pbmcapply::pbmclapply,cores = 12)->read_replicates_output_no_scn_raw)
# 

## ----select_similar_samples,eval=F--------------------------------------------
# read_replicates_output_no_scn_raw$sample_mat %>% as.matrix() %>% Rfast::colsums() ->rs_no_scn_raw
# lapply(1:length(unique(conditions)),function(i) {which(conditions==unique(conditions)[i])[which(rank(rs_no_scn_raw[which(conditions==unique(conditions)[i])])>3)]}) %>% setNames(unique(conditions)) ->top_samples
# lapply(1:length(unique(conditions)),function(i) {which(conditions==unique(conditions)[i])[which(rank(rs_no_scn_raw[which(conditions==unique(conditions)[i])])<=3)]}) %>% setNames(unique(conditions)) ->bottom_samples
# 
# c(bottom_samples$trt,top_samples$ctrl)->selected_samples
# read_replicates_output_no_scn_selected=read_replicates_output_no_scn_raw
# read_replicates_output_no_scn_selected$sample_mat[,..selected_samples]->read_replicates_output_no_scn_selected$sample_mat
# 

## ----equalize_reads,eval=F----------------------------------------------------
# downsample_rr(read_replicates_output_no_scn_selected)->read_replicates_output_no_scn_ds
# read_replicates_output_no_scn_raw$conditions[selected_samples]->read_replicates_output_no_scn_ds$conditions
# #notice how all of the samples have exactly the same read counts at this point.
# 

## ----filter,eval=F------------------------------------------------------------
# basepairC::normalizeFilterReadCounts(sample_mat = read_replicates_output_no_scn_ds$sample_mat,conditions = read_replicates_output_no_scn_ds$conditions,loc1 = read_replicates_output_no_scn_ds$loc1,loc2 = read_replicates_output_no_scn_ds$loc2,nrow_output = read_replicates_output_no_scn_ds$rep_dim[1],ncol_output = read_replicates_output_no_scn_ds$rep_dim[2],output_type = "corrected_counts",filter = T,underscored_positions_col = read_replicates_output_no_scn_ds$underscored_positions_col,underscored_positions_row = read_replicates_output_no_scn_ds$underscored_positions_row,norm_factor_type = "DESeq2",filter_type = "diff")->normalizeFilterReadCounts_output_counts_no_scn_deseq2_diff_ds

## ----save_and_load_normalized_matrices,eval=F,echo=F--------------------------
# #saveRDS(normalizeFilterReadCounts_output_counts_no_scn_deseq2_diff_ds,file = "normalizeFilterReadCounts_output_counts_no_scn_deseq2_diff_ds.rds")
# #normalizeFilterReadCounts_output_counts_no_scn_deseq2_diff_ds=readRDS(file = "normalizeFilterReadCounts_output_counts_no_scn_deseq2_diff_ds.rds")

## ----get_matrices,eval=F------------------------------------------------------
# readRDS(system.file(package="basepairC","data/mnase_bias_6bp.rds"))[1:6000,1:6000]->mnase_bias_matrix
# 
# normalizeFilterReadCounts_output_counts_no_scn_deseq2_diff_ds$sum_mat_list$ctrl[1:6000,1:6000] ->ctrl_mat_no_scn_deseq2_diff_ds
# normalizeFilterReadCounts_output_counts_no_scn_deseq2_diff_ds$sum_mat_list$trt[1:6000,1:6000] ->trt_mat_no_scn_deseq2_diff_ds

