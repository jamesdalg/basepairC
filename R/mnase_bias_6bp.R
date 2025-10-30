#' Calculate MNase bias in a Sliding Window
#'
#' This function calculates the MNase bias within a sliding window across specified genomic positions on a chromosome.
#'
#' @param chromosome The chromosome for which to calculate MNase bias.
#' @param position A vector of positions on the chromosome to calculate MNase bias for.
#' @param genome The genome object from which to retrieve sequence data. By default, this is set to the mouse genome (BSgenome.Mmusculus.UCSC.mm10).
#' @param window_size The size of the sliding window for GC content calculation. Default is 5.
#' @param n_cores The number of cores to use for parallel processing. Default is half of the available cores detected on the system.
#'
#' @return A vector of GC content values for each position, calculated within the specified sliding window.
#'
#' @examples
#' # Calculate GC content for positions 100 to 200 on chromosome 1
#' mnase_bias_test <- mnase_bias_win(
#'   chromosome = "chr1",
#'   position = 100:200,
#'   genome = BSgenome.Mmusculus.UCSC.mm10,
#'   window_size = 6
#' )
#'
#' @importFrom Biostrings getSeq
#' @importFrom data.table data.table
#' @importFrom parallel detectCores
#' @importFrom GenomicRanges GRanges
#' @importFrom pbmcapply pbmclapply
#' @export
mnase_bias_6bp <- function(chromosome, position, genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10, window_size = 6, n_cores = parallel::detectCores() / 2,mnase_model=system.file(package="basepairC","extdata/mnase_model_6bp.rds")) {
  #browser()
  if(class(mnase_model) %in% "character"){
    if(stringr::str_detect(string = tolower(mnase_model),pattern = ".rds")){
      mnase_model <- readRDS(mnase_model)
      mnase_model[,c(1,4)]->mnase_model
      colnames(mnase_model)=c("seq","bias")
      mnase_vec=mnase_model$bias
      names(mnase_vec)=mnase_model$seq
    } else{
    if(stringr::str_detect(string = tolower(mnase_model),pattern = ".rda")){
      load(mnase_model)
      mnase_model_6bp[,c(1,4)]->mnase_model
      colnames(mnase_model)=c("seq","bias")
      mnase_vec=mnase_model$bias
      names(mnase_vec)=mnase_model$seq
    }
    }
    
  }
  # if(is.null(mnase_model)){
  #   mnase_model_6bp[,c(1,4)]->mnase_model
  #   colnames(mnase_model)=c("seq","bias")
  #   mnase_vec=mnase_model$bias
  #   names(mnase_vec)=mnase_model$seq
  # }
  data.table::data.table(chr = unique(chromosome), start = (position-2), end = (position+3)) -> gene_locus_dt_bp
  #data.table::data.table(chr = unique(chromosome), start = (position-2), end = (position+3)) %>% GRanges() -> gene_locus_gr
  #Biostrings::getSeq(genome, gene_locus_gr+window_size) -> gene_locus_seq
  #data.table::data.table(chr = unique(as.character(chromosome)), start = min(position)-2, end=max(position)+3) %>%
#    dplyr::mutate(seq = gene_locus_seq %>% as.character() %>% strsplit("") %>% unlist(), end = start) %>%
    #dplyr::select(chr, start, end) -> gene_locus_dt_bp
   #browser()
  pbmcapply::pbmclapply(1:(nrow(gene_locus_dt_bp)), function(i) {
    # print(i)
    #min_window_pos <- max(i-floor(window_size / 2-0.5), 1-floor(window_size / 2-0.5))#, #i-floor(window_size / 2+0.5)+(window_size %% 2)+1 #i-floor(window_size / 2+0.5)+(window_size %% 2)
    #max_window_pos <- min(i + floor(window_size / 2+0.5)-(window_size %% 2), (nrow(gene_locus_dt_bp)-window_size)+floor(window_size/2)-window_size%%2 )
    gene_locus_dt_bp[i,] %>% GRanges() %>% Biostrings::getSeq(genome, .) %>% as.character()-> window_seq
    mnase_vec[window_seq] -> mnase_bias_win
    if(is.na(mnase_vec[window_seq])){
      mnase_bias_win=0
    }
    #sum(window_seq %in% c("G", "C")) / length(window_seq) -> gc_content
    
    #  if(any(is.na(mnase_bias_win))){
    #   which(is.na(mnase_bias_win))->na_pos
    #   mnase_bias_win[na_pos]=0
    # }
    return(mnase_bias_win)
  }, mc.cores = n_cores) %>% unlist() -> mnase_bias_vec
  return(mnase_bias_vec)
}
