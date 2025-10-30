#' Get CTCF Segments
#'
#' This function retrieves CTCF binding sites from a genomic region, segments a sample matrix list,
#' and returns the corresponding breakpoints.
#'
#' @param chr A character string indicating the chromosome (default is "chr15").
#' @param start An integer representing the start position of the genomic region (default is 61982445).
#' @param end An integer representing the end position of the genomic region (default is 61988445).
#' @param genome_name A character string indicating the genome version (default is "mm10").
#' @param sample_mat_list A list of matrices representing sample data.
#'
#' @return A list containing breakpoint positions for rows and columns, and their transposed versions.
#' @importFrom AnnotationHub AnnotationHub subset
#' @importFrom GenomicRanges GRanges reduce
#' @importFrom IRanges IRanges  subsetByOverlaps
#' @importFrom dplyr bind_rows mutate pull
#' @importFrom pbapply pblapply
#' @importFrom utils data
#' @import CTCF
#' @export
get_ctcf_segments=function(chr="chr15",start=61982445,end=61988445,genome_name="mm10",sample_mat_list=NULL,return_gr=F,overlap_gap=0){
  #browser()
  if(is.null(sample_mat_list)){stop("sample_mat_list is null")}
  row_underscored_pos=sample_mat_list[[1]] %>% rownames()
  col_underscored_pos=sample_mat_list[[1]] %>% colnames()
  ah <- AnnotationHub::AnnotationHub()
  subset(ah, ah$genome == genome_name & ah$preparerclass=="CTCF")->genome_ah
  pbapply::pblapply(1:length(genome_ah),function(i){subsetByOverlaps(genome_ah[[i]], GRanges(chr, IRanges(start, end))) %>% as.data.frame() %>%  dplyr::mutate(source=genome_ah$title[i])})->ctcf_sites_l
  do.call(dplyr::bind_rows,ctcf_sites_l)->ctcf_sites
  (ctcf_sites %>% GenomicRanges::GRanges()+overlap_gap) %>% GenomicRanges::reduce() ->reduced_gr
  reduced_gr %>% as.data.frame() %>% dplyr::mutate(rng=paste0(seqnames,"_",start,"_",start)) %>% dplyr::pull(rng) -> ctcf_seg_underscored_ranges_start
  reduced_gr %>% as.data.frame() %>% dplyr::mutate(rng=paste0(seqnames,"_",end,"_",end)) %>% dplyr::pull(rng) -> ctcf_seg_underscored_ranges_end
  combined_positions=unique(c(ctcf_seg_underscored_ranges_start,ctcf_seg_underscored_ranges_end))
  which(col_underscored_pos %in% combined_positions)->ctcf_seg_num_col
  which(row_underscored_pos %in% combined_positions)->ctcf_seg_num_row
  seg_object=list(breakpoints_col=ctcf_seg_num_col,breakpoints_row=ctcf_seg_num_row,t_breakpoints_col=ctcf_seg_num_row,t_breakpoints_row=ctcf_seg_num_col)
  if(return_gr){return(list(segments=seg_object,gr=reduced_gr))  } else {return(seg_object)}

}