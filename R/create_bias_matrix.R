#' Create a MNase Bias Matrix
#'
#' This function generates a bias matrix based on the MNase-seq bias profile for a given genomic region. 
#'
#' @param chr A character string specifying the chromosome of the region of interest. Default is `"chr15"`.
#' @param start An integer specifying the start position of the region of interest. Default is `61982445`.
#' @param end An integer specifying the end position of the region of interest. Default is `61988445`.
#' @param genome A BSgenome object representing the genome of interest. Default is `BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10`.
#' @param map_fn A character string specifying the file path to the map file or `"mm10_k24"` to download it automatically. Default is the `k24.umap.bedgraph.gz` file from the `basepairC` package.
#' @param nrow An integer specifying the number of rows in the output matrix. Default is `6000`.
#' @param ncol An integer specifying the number of columns in the output matrix. Default is `6000`.
#'
#' @return A numeric matrix of size `nrow` x `ncol`, where each element represents the MNase-seq bias for a given genomic position.
#'
#' @details 
#' Using the exact genomic sequence of the specified region,
#' constructs a bias matrix by multiplying a bias row vector with its transpose,
#'  resulting in a symmetric matrix representing the MNase-seq bias.
#'
#' @examples
#' \dontrun{
#' bias_matrix <- create_bias_matrix(chr="chr15", start=61982445, end=61988445)
#' }
#'
#' @export
create_bias_matrix=function(chr="chr15",start=61982445,end=61988445,genome=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10,map_fn=system.file(package="basepairC","data/mm10/k24.umap.bedgraph.gz"),nrow=6000,ncol=6000){
  if(map_fn=="mm10_k24"){
    message_parallel("Downloading the map file")
    if(
      (file.exists("mm10_k24.umap.bedgraph.gz"))
    ){
      #unzip the mm10_k24.umap.bedgraph.gz
      R.utils::gunzip("mm10_k24.umap.bedgraph.gz",remove=F)->map_fn
    } else {
      download.file("https://bismap.hoffmanlab.org/raw/mm10/k24.umap.bedgraph.gz","mm10_k24.umap.bedgraph.gz")
      R.utils::gunzip("mm10_k24.umap.bedgraph.gz",remove=F)->map_fn
    }
  }
  if(stringr::str_detect(string = map_fn,pattern=".gz")){
    if(!file.exists(gsub(map_fn,".gz",""))){
      message_parallel("Unzipping the map file")
      R.utils::gunzip(map_fn,remove=F,temporary=T,overwrite=T)->map_fn}
    #map_fn=gsub(".gz","",map_fn)
  }
  message_parallel("Creating the bias dt")
  create_dt_with_map_and_gc(chr_val=chr,start_val=start,end_val=end,genome_val=genome,map_fn=map_fn)->dt_with_map_and_gc
  message_parallel("Creating the bias matrix")
  matrix(rep(dt_with_map_and_gc$mnase_bias_6bp,length(dt_with_map_and_gc$mnase_bias_6bp)),nrow=length(dt_with_map_and_gc$mnase_bias_6bp),byrow=T) ->mnase_bias_row
  mnase_bias_row %>% t()->mnase_bias_col
  mnase_bias_matrix=(mnase_bias_row*mnase_bias_col)
  rownames(mnase_bias_matrix)<-colnames(mnase_bias_matrix)<-paste0(dt_with_map_and_gc$chr,"_",dt_with_map_and_gc$start,"_",dt_with_map_and_gc$end)
  mnase_bias_matrix[1:nrow,1:ncol]->mnase_bias_matrix
  file.remove(map_fn)
  return(mnase_bias_matrix)
}
