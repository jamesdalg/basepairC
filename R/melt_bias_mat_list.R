#' Melt a List of Data Tables Containing Genomic Mapping and GC Data
#'
#' This function takes a list of data tables, each containing genomic mapping and GC data, and melts them into a long format. It optionally uses parallel processing to speed up the operation. The first element of the list retains the genomic coordinates, while the rest are concatenated without the coordinates.
#'
#' @param dt_with_map_and_gc_list A list of data tables, where each element contains genomic mapping and GC data. The names of the list elements are used as column names in the melted data table.
#' @param ncores Integer. The number of cores to use for parallel processing. Default is 1.
#'
#' @return A melted data table containing genomic mapping and GC data in long format. The first table retains the genomic coordinates (`chr1`, `start1`, `end1`, `chr2`, `start2`, `end2`), while the remaining tables are concatenated by their values.
#'
#' @details
#' - The function uses `pbapply::pblapply` for parallel processing of each data table in the list.
#' - The first data table in the list is split into columns `chr1`, `start1`, `end1`, `chr2`, `start2`, `end2` using `tstrsplit`, while subsequent tables only include the melted values without the coordinates.
#' - The column names `loc1` and `loc2` are used to represent the genomic locations in the first data table, which are split by the separator `_` or space to extract the chromosome and start/end positions.
#' 
#' @import data.table pbapply reshape2 dplyr
#' 
#' @examples
#' # Example usage:
#' dt_list <- list(dt1 = data.table(matrix(rnorm(25), nrow=5)),
#'                 dt2 = data.table(matrix(rnorm(25), nrow=5)))
#' names(dt_list) <- c("map_gc_1", "map_gc_2")
#' result <- melt_bias_mat_list(dt_list, ncores = 2)
#' head(result)
#' 
#' @export
melt_bias_mat_list=function(dt_with_map_and_gc_list,ncores=1) {
  data.table::setDTthreads(ncores)
  pbapply::pblapply(1:length(dt_with_map_and_gc_list), function(i) {
    reshape2::melt(dt_with_map_and_gc_list[[i]]) %>% data.table::as.data.table()->melted_dt
    colnames(melted_dt)=c("loc1","loc2",names(dt_with_map_and_gc_list)[i])
    if(i!=1){melted_dt %>% dplyr::select(-loc1,-loc2)->melted_dt} else {
      melted_dt[, c("chr1", "start1", "end1") := tstrsplit(loc1, split = "_| ")]
      melted_dt[, c("chr2", "start2", "end2") := tstrsplit(loc2, split = "_| ")]
      melted_dt[, `:=`(start1 = as.numeric(start1), start2 = as.numeric(start2))]
    }
    return(data.table::as.data.table(melted_dt))
  })   %>% do.call(cbind,.)  ->dt_with_map_and_gc_list_melted
}
