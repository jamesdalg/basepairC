library(magrittr)
setwd("/code/R/allele_spec_comparisons/")
skew_files=list.files("./data/alelleic_combined/rs2793109_matrices/raw",pattern=glob2rx("*sk*.txt"),full.names=TRUE,recursive=T)
dir_df=data.frame(
  skew_file=skew_files,
  dir=dirname(skew_files),
  dataset=basename(dirname(skew_files))
)
d=1
pbapply::pblapply(1:rev(1:nrow(dir_df)),function(d){ 
  tryCatch({
  basepairC.core::message_parallel(paste0("Processing dataset ",dir_df$dataset[d]))
small_matrix_fns=setdiff(list.files(dir_df$dir[d],full.names=TRUE,pattern=glob2rx("*.txt")),skew_files)

sample_sheet=data.frame(
  sample_name=gsub(".txt","",basename(small_matrix_fns)),
  small_matrix_fn_base=basename(small_matrix_fns),
  small_matrix_fn_full=small_matrix_fns
) |> tidyr::separate(sample_name,into=c("rsid", "allele","P1","P2","cell_type","sample","sample2"),sep="_|-",remove = F) %>% tidyr::unite("sample",sample:sample2,sep="|") |> dplyr::select(small_matrix_fn_full,allele,sample,rsid) %>% dplyr::mutate(rsid_allele=paste0(rsid,"_",allele)) 
unique_rsids=unique(sample_sheet$rsid)
unique_rsids=setdiff(unique_rsids,"rs1474379684")
i=1
i=which(unique_rsids=="rs11187157")
i=which(unique_rsids=="rs1474379684")
i=which(unique_rsids=="rs11187157")
data.table::setDTthreads(1)
i=1
j=1





options(mc.cores=parallel::detectCores()/16) 
pbmcapply::pbmclapply(1:length(unique_rsids),function(i){
  basepairC.core::message_parallel(paste0("Processing rsid ",unique_rsids[i]))
  tryCatch({
data.table::setDTthreads(1)
basepairC.core::read_replicates(replicate_files = sample_sheet %>% dplyr::filter(rsid==unique_rsids[i]) %>% dplyr::pull(small_matrix_fn_full),
                                conditions = sample_sheet %>% dplyr::filter(rsid==unique_rsids[i]) %>% dplyr::pull(allele),
                                thresh=0,
                                cap=Inf,
                                scn = F,
                                python = T,
                                apply_fun=pbmcapply::pbmclapply,
                                cores = 1,
                                reorder=F )->read_replicates_output_no_scn_raw
zero_inf_mat=read_replicates_output_no_scn_raw$sample_mat %>% as.matrix()
zero_inf_mat[zero_inf_mat>0]=1
Rfast::colsums(zero_inf_mat)->non_zero_counts_per_col
nonzero_percent_per_col=non_zero_counts_per_col/nrow(zero_inf_mat)

read_replicates_output_no_scn_raw %>% basepairC.core::cpm_filter_rc_norm_dgel(filter_type="none",norm_factor_type = "TMM")->mat_list_no_scn_raw_dgel
mat_list_no_scn_raw_dgel$design=model.matrix(~1+as.factor(sample_sheet$allele))


mat_list_no_scn_raw_dgel$counts->count_mat
colnames(count_mat)=paste0(unique_rsids[i],
                                                 "_",
                                                 sample_sheet %>% dplyr::filter(rsid==unique_rsids[i]) %>% dplyr::pull(allele),
                                                 "_",
                                                 sample_sheet %>% dplyr::filter(rsid==unique_rsids[i]) %>% dplyr::pull(sample))
cbind(count_mat,mat_list_no_scn_raw_dgel$loc1,mat_list_no_scn_raw_dgel$loc2) %>% data.table::as.data.table() %>% data.table::melt.data.table(id.vars=c('start1','end1','chr1','start2','end2','chr2'),
                                                                                                                                             variable.name = "rsid_allele_sample",
                                                                                                                                             value.name = "count") -> mat_long_df
mat_long_df$rsid_allele_sample=as.character(mat_long_df$rsid_allele_sample)
data.table::tstrsplit(mat_long_df$rsid_allele_sample, "_", keep = 1:3, names = c("rsid_allele", "overall_sample_num","sample")) %>% as.data.frame()-> rsid_sample_allele_split_df
colnames(rsid_sample_allele_split_df)=c("rsid","allele","sample")
mat_long_df_with_metadata=cbind(rsid_sample_allele_split_df,mat_long_df)
skew_vars_df=data.table::fread(dir_df$skew_file[d],col.names=c("rsid","chr","pos","ref","alt"))
j=1
pbapply::pblapply(1:10,function(j){ 
  tryCatch({

data.table::merge.data.table(mat_long_df_with_metadata,skew_vars_df,by="rsid")->mat_long_merged_df
mat_long_merged_df$alt_binary=ifelse(mat_long_merged_df$allele==mat_long_merged_df$alt,1,0)
mat_long_merged_df %>% dplyr::select(alt_binary,sample) %>% dplyr::distinct()-> true_trt_df
mat_long_merged_df = as.data.frame(mat_long_merged_df)
alt_sample_df=mat_long_merged_df %>% dplyr::select(alt_binary,sample) %>% dplyr::distinct()
set.seed(j)

rand_trt=sample(c(rep(1, round(nrow(alt_sample_df)/4)), rep(0, round(nrow(alt_sample_df)/4))),size = round(nrow(alt_sample_df)/2))
set.seed(j*1000)
rand_control=sample(c(rep(1, round(nrow(alt_sample_df)/4)), rep(0, round(nrow(alt_sample_df)/4))),size = round(nrow(alt_sample_df)/2))
sham_trt_df=data.frame(
  sample=unique(mat_long_merged_df$sample),
  sham_trt=c(rand_trt,rand_control)
)
if(all(dim(sham_trt_df)==dim(true_trt_df)) ){
   basepairC.core::message_parallel(paste0("sham_trt_df and true_trt_df have matching dimensions for rsid ",unique_rsids[i],", seed ",j)) } else {
   basepairC.core::message_parallel(paste0("sham_trt_df and true_trt_df DO NOT have matching dimensions for rsid ",unique_rsids[i],", seed ",j)) 
  return(
   data.frame()
  )
}
mat_long_merged_df_sham_merged=dplyr::left_join(mat_long_merged_df,sham_trt_df,by="sample",relationship="many-to-many")
mat_long_merged_df_sham_merged$sham_trt_factor=ifelse(mat_long_merged_df_sham_merged$sham_trt==1,"sham_trt","sham_ctrl")
mat_long_merged_df_sham_merged %>% dplyr::filter(mat_long_merged_df_sham_merged$rsid==unique_rsids[i],alt_binary==0)->ref_data_df
mat_long_merged_df_sham_merged %>% dplyr::filter(mat_long_merged_df_sham_merged$rsid==unique_rsids[i],alt_binary==1)->alt_data_df
mat_long_merged_df_sham_merged %>% dplyr::filter(mat_long_merged_df_sham_merged$rsid==unique_rsids[i])->full_data_df
full_data_df %>% dplyr::mutate(alt_factor=factor(ifelse(alt_binary==1,"alt","ref"),levels=c("ref","alt"))) ->full_data_df
full_data_df$alt_factor=relevel(full_data_df$alt_factor,ref="ref")
starting_values=NULL
tryCatch({
    starting_values <- basepairC.core::bootstrap_gam_parallel_speedglm(
      data = ref_data_df,
      formula = count ~  sham_trt_factor, 
      family = quasipoisson(),
      fraction = 0.01, 
      n_bootstraps = 5
    ) %>% .$avg_coefs 
  }, error = function(e) {
    basepairC.core::message_parallel(paste0("Error in bootstrap for rsid ",unique_rsids[i], ": ", e))
    starting_values <- NULL
  })
speedglm::speedglm(count ~  sham_trt_factor,trace=T, data=ref_data_df,family=quasipoisson(),start=starting_values)->qp_model_sham_ref
starting_values=NULL
tryCatch({
    starting_values <- basepairC.core::bootstrap_gam_parallel_speedglm(
      data = alt_data_df,
      formula = count ~  sham_trt_factor, 
      family = quasipoisson(),
      fraction = 0.01, 
      n_bootstraps = 5
    ) %>% .$avg_coefs 
  }, error = function(e) {
    basepairC.core::message_parallel(paste0("Error in bootstrap for rsid ",unique_rsids[i], ": ", e))
    starting_values <- NULL
  })

speedglm::speedglm(count ~  sham_trt_factor,trace=T, data=alt_data_df,family=quasipoisson(),start=starting_values)->qp_model_sham_alt
broom.mixed::tidy(qp_model_sham_ref,parametric=T)->qp_model_tidy_ref
qp_model_tidy_ref$group="ref_sham"
broom.mixed::tidy(qp_model_sham_alt,parametric=T)->qp_model_tidy_alt
qp_model_tidy_alt$group="alt_sham"

qp_model_tidy=rbind(qp_model_tidy_ref,qp_model_tidy_alt)
qp_model_tidy$random_seed=j
return(qp_model_tidy)
},error=function(e){
  basepairC.core::message_parallel(paste0("Error with rsid ",unique_rsids[i]," replicate ",j,": ",e))
  return( data.frame() )
}) }) %>% do.call(rbind,.)->qp_model_tidy_replicates
qp_model_tidy_replicates %>% dplyr::group_by(term) %>% dplyr::summarise(est_abs_mean=mean(abs(estimate),na.rm=T),
                                                                        est_abs_sd=sd(abs(estimate),na.rm=T),
                                                                        est_mean=mean(estimate,na.rm=T),
                                                                        mean_se=mean(std.error,na.rm=T), 
                                                                        mean_abs_stat=mean(abs(statistic),na.rm=T),
                                                                        mean_p=mean(p.value,na.rm=T),
                                                                        est_sd=sd(estimate,na.rm=T)) %>%
  dplyr::mutate(logp_mean_stat=Rmpfr::pnorm(mean_abs_stat,lower.tail = F,log.p=T))->qp_model_tidy_summary
data.table::merge.data.table(mat_long_df_with_metadata,skew_vars_df,by="rsid")->mat_long_merged_df_full

mat_long_merged_df_full$alt_binary=ifelse(mat_long_merged_df_full$allele==mat_long_merged_df_full$alt,1,0)
mat_long_merged_df_full = as.data.frame(mat_long_merged_df_full)
mat_long_merged_df_full %>% dplyr::mutate(alt_factor=factor(ifelse(alt_binary==1,"alt","ref"),levels=c("ref","alt"))) ->full_data_df_true
full_data_df_true$alt_factor=relevel(full_data_df_true$alt_factor,ref="ref")   

starting_values=NULL
tryCatch({
    starting_values <- basepairC.core::bootstrap_gam_parallel_speedglm(
      data = full_data_df_true,
      formula = count ~ alt_factor, 
      family = quasipoisson(),
      fraction = 0.01, 
      n_bootstraps = 5
    ) %>% .$avg_coefs 
  }, error = function(e) {
    basepairC.core::message_parallel(paste0("Error in bootstrap for rsid ",unique_rsids[i], ": ", e))
    starting_values <- NULL
  })
    
speedglm::speedglm(count ~ alt_factor,trace=T, data=full_data_df_true,family=quasipoisson(),start=starting_values)->qp_model_true_speedglm
broom.mixed::tidy(qp_model_true_speedglm,parametric=T)->qp_model_tidy_true
qp_model_tidy_summary$analysis="sham"
qp_model_tidy_true$analysis="true"
qp_model_tidy_complete=plyr::rbind.fill(qp_model_tidy_true,qp_model_tidy_summary)
qp_model_tidy_complete$rsid=unique_rsids[i]
qp_model_tidy_complete$nonzero_percent_per_col=paste0(round(nonzero_percent_per_col*100,digits = 6),collapse=",")
qp_model_tidy_complete$nonzero_counts_per_col=paste0(non_zero_counts_per_col,collapse=",")
qp_model_tidy_complete$failed=F
return(qp_model_tidy_complete)
  },error=function(e){
    basepairC.core::message_parallel(paste0("Error with rsid ",unique_rsids[i],": ",e))
    return(
      data.frame(effect=NA,component=NA,group=NA,term=NA,estimate=NA,std.error=NA,statistic=NA,p.value=NA,rsid=unique_rsids[i],fit_message=NA,failed=TRUE)
      )
  })
}) %>% data.table::rbindlist(fill=T)->model_df_complete_glmmtmb_ziqp #

dir.create("/code/R/allele_spec_comparisons/results/",showWarnings = FALSE)
dir.create("/code/R/allele_spec_comparisons/results/intermediate_all_jan_17/",showWarnings = FALSE)
dir.create("/code/R/allele_spec_comparisons/results/merged_all_jan_17_2026/",showWarnings = FALSE)
dir.create("/code/R/allele_spec_comparisons/dataset_errors_jan_17/",showWarnings = FALSE)
data.table::fwrite(x = model_df_complete_glmmtmb_ziqp,paste0("/code/R/allele_spec_comparisons/results/intermediate_all_jan_17/allele_spec_comparison_neg_control_speedglm_qp_results_glm__intercept_with_sham_",dir_df$dataset[d],".tsv"),sep="\t")
model_df_complete_glmmtmb_ziqp %>% dplyr::select(-est_abs_mean,-est_abs_sd,-est_mean,-est_sd,-mean_se,-mean_abs_stat,-mean_p,-logp_mean_stat) %>% dplyr::filter(analysis=="true") %>% data.table::fwrite(paste0("/code/R/allele_spec_comparisons/results/intermediate_all_jan_17/allele_spec_comparison_neg_control_speedglm_qp_results_glm__intercept_true_analysis_",dir_df$dataset[d],".tsv"),sep="\t")
model_df_complete_glmmtmb_ziqp %>% dplyr::select(-estimate,-std.error,-statistic,-p.value) %>% dplyr::filter(analysis=="sham") %>% data.table::fwrite(paste0("/code/R/allele_spec_comparisons/results/intermediate_all_jan_17/allele_spec_comparison_neg_control_speedglm_qp_results_glm__intercept_sham_analysis_",dir_df$dataset[d],".tsv"),sep="\t")

model_df_complete_glmmtmb_ziqp %>% dplyr::select(-estimate,-std.error,-statistic,-p.value) %>% dplyr::filter(analysis=="sham")   %>% dplyr::filter(term=="sham_trt_factorsham_trt") %>% dplyr::mutate(analysis=as.character(analysis)) %>% dplyr::select(term,est_abs_mean,est_mean,est_abs_sd,rsid) %>% dplyr::inner_join(
  model_df_complete_glmmtmb_ziqp %>% dplyr::select(-est_abs_mean,-est_abs_sd,-est_mean,-est_sd,-mean_se,-mean_abs_stat,-mean_p,-logp_mean_stat) %>% dplyr::filter(analysis=="true")     %>% dplyr::filter(term=="alt_factoralt") %>% dplyr::mutate(analysis=as.character(analysis)) %>% dplyr::select(rsid,term,estimate,std.error,p.value,nonzero_percent_per_col,nonzero_counts_per_col),
  by="rsid"
) -> sham_trt_merged
data.table::fwrite(sham_trt_merged,paste0("/code/R/allele_spec_comparisons/results/merged_all_jan_17_2026/allele_spec_comparison_neg_control_speedglm_qp_results_glm_no_sample_intercept_sham_trt_vs_true_alt_effect_",dir_df$dataset[d],"_NO_SAMPLE.tsv"),sep="\t")
return(sham_trt_merged)
  },error=function(e){
    basepairC.core::message_parallel(paste0("Error with dataset ",dir_df$dataset[d],": ",e))
    writeLines(text = paste0(e),con=paste0("/code/R/allele_spec_comparisons/dataset_errors_jan_17/",dir_df$dataset[d],".error"))
    return(data.frame())
  })
}) %>% data.table::rbindlist(fill=TRUE)->all_datasets_sham_trt_vs_true_alt_effect
data.table::fwrite(all_datasets_sham_trt_vs_true_alt_effect,"/code/R/allele_spec_comparisons/results/merged_all_jan_17_2026/allele_spec_comparison_neg_control_speedglm_qp_results_glm_intercept_sham_trt_vs_true_alt_effect_ALL_datasets_NO_SAMPLE.tsv",sep="\t")

