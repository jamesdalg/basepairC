#' Convert Cell IDs to Genomic Interactions
#'
#' This function takes a vector of cell IDs and converts them into a 
#' \code{GenomicInteractions} object. The cell IDs are expected to be 
#' formatted as strings with chromosome and genomic coordinates separated 
#' by underscores or spaces. The function splits the cell IDs into 
#' components, creates two \code{GRanges} objects (one for each genomic 
#' anchor), and combines them into a \code{GenomicInteractions} object.
#'
#' @param cell_ids A character vector of cell IDs, where each element is 
#' expected to contain chromosome and coordinate information separated by 
#' underscores or spaces. The expected format for each ID is: 
#' \code{chr1_start1_end1_chr2_start2_end2}.
#'
#' @return A \code{GenomicInteractions} object representing the genomic 
#' interactions defined by the input cell IDs.
#'
#' @importFrom reshape2 colsplit
#' @importFrom dplyr rename
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicInteractions GenomicInteractions
#'
#' @examples
#' cell_ids <- c("chr1_100_200_chr2_300_400", "chr1_500_600_chr2_700_800")
#' interactions <- cell_ids_to_gint(cell_ids)
#'
#' @export
cell_ids_to_gint=function(cell_ids){
  cell_ids %>% reshape2::colsplit(pattern="_| ",names=c("chr1","start1","end1","chr2","start2","end2")) ->cell_id_df
  GenomicInteractions::GenomicInteractions(anchor1 = cell_id_df[,1:3] %>% dplyr::rename(chr=chr1,start=start1,end=end1) %>% GenomicRanges::GRanges(),anchor2 = cell_id_df[,4:6] %>% dplyr::rename(chr=chr2,start=start2,end=end2) %>% GenomicRanges::GRanges())
}
gint_to_cell_ids<-function(gint){
  gint %>% as.data.frame() %>% dplyr::mutate(cell_id=paste0(seqnames1,"_",start1,"_",end1,"_",seqnames2,"_",start1,"_",end1)) %>% dplyr::pull(cell_id)
}
mirror_cells=function(df,cell_id_name="cell_id"){
  df %>% tidyr::separate(cell_id_name,sep="_| ",into=c("cell_chr1","cell_start1","cell_end1","cell_chr2","cell_start2","cell_end2")) %>% dplyr::mutate(cell_id=paste0(cell_chr1,"_",cell_start1,"_",cell_end1," ",cell_chr2,"_",cell_start2,"_",cell_end2)) ->split_df
  split_df %>% dplyr::mutate(cell_id=paste0(cell_chr1,"_",cell_start1,"_",cell_end1," ",cell_chr2,"_",cell_start2,"_",cell_end2)) %>% dplyr::select(-cell_start1,-cell_start2,-cell_chr1,-cell_chr2,-cell_end1,-cell_end2) ->split_df_orig
  split_df %>% dplyr::rename(new_cell_chr1=cell_chr2,new_cell_start1=cell_start2,new_cell_end1=cell_end2,new_cell_chr2=cell_chr1,new_cell_start2=cell_start1,new_cell_end2=cell_end1) %>% 
    dplyr::mutate(cell_id=paste0(new_cell_chr1,"_",new_cell_start1,"_",new_cell_end1," ",new_cell_chr2,"_",new_cell_start2,"_",new_cell_end2)) %>% dplyr::select(-new_cell_chr1,-new_cell_chr2,-new_cell_start1,-new_cell_start2,-new_cell_end1,-new_cell_end2)->split_df_mirror
  return(dplyr::bind_rows(split_df_orig,split_df_mirror) %>% dplyr::distinct())
}