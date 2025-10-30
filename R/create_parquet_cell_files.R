create_parquet_cell_files <- function(parquet_dataset_filename,parquet_dir){
  system.time({arrow::read_parquet(parquet_dataset_filename)->in_memory_parquet_dataset})
  setkey(in_memory_parquet_dataset, lasso_cell_id)
  #create ../dev/parquet_cell directory.
  dir.create(parquet_dir,showWarnings = F,recursive = T)
  #pull out the cells one at a time and write to a parquet file
  in_memory_parquet_dataset %>% dplyr::select(lasso_cell_id) %>% dplyr::distinct() %>% dplyr::collect() %>% unlist() %>% unique()->cell_ids
  pbapply::pblapply(1:length(cell_ids), function(i) {
    #get the cell data
    message_parallel(paste0("writing cell:",cell_ids[i]))
    tryCatch({
      in_memory_parquet_dataset %>% dplyr::filter(lasso_cell_id==cell_ids[i]) %>% dplyr::collect()->cell_df
      arrow::write_parquet(cell_df, paste0(parquet_dir,cell_ids[i],".parquet"))
      return(data.frame())
    },error=function(e){warning(paste0(cell_ids[i]," failed to write"));
      message_parallel(paste0(cell_ids[i]," failed to write",i));
      return(data.frame())
      ;})
  })
}