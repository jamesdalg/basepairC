#' Calculate MNase bias in a Sliding Window
#'
#' This function calculates the MNase bias within a sliding window across specified genomic positions on a chromosome.
#'
#' @param chromosome The chromosome for which to calculate MNase bias.
#' @param position A vector of positions on the chromosome to calculate MNase bias for.
#' @param genome The genome object from which to retrieve sequence data. By default, this is set to the mouse genome (BSgenome.Mmusculus.UCSC.mm10).
#' @param kmer_size The size of the k-mer to use for calculating MNase bias. Default is 6.
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
mnase_bias_kbp <- function(chromosome, position, genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10, kmer_size = 6, n_cores = parallel::detectCores() / 2,mnase_model=system.file(package="basepairC","data/mnase_model_6bp.rds"),return_seq_model=F,bias_type="ratio") {
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
      #colnames(standard_6bp_model)=c("seq","obs","exp","o_e_ratio")
    }
      #browser()
      #
      data.table::fread(mnase_model,skip=1)->mnase_model #
      colnames(mnase_model)=c("index", "seq","plus_exp","minus_exp","plus_obs","minus_obs")  #exp is from the fasta, obs is from the bam.
      mnase_model %>% dplyr::select(-index) %>% dplyr::mutate(o_e_ratio=ifelse((plus_exp+minus_exp)>0,
                                              (plus_obs+minus_obs)/(plus_exp+minus_exp),0),
                                              relative_bias=ifelse((plus_exp+minus_exp)>0,
                                            ((plus_obs+minus_obs)-(plus_exp+minus_exp))/(plus_exp+minus_exp),0))->mnase_model
      if(bias_type=="ratio"){
        mnase_vec=mnase_model$o_e_ratio
      } 
      if(bias_type=="relative"){
        mnase_vec=mnase_model$relative_bias
      }
      
      names(mnase_vec)=mnase_model$seq
      if(return_seq_model){
        return(mnase_model)
      }
      #please note, this is a relative bias measure and not the traditional bias=o-e.
    }
    
  }
  # if(is.null(mnase_model)){
  #   mnase_model_6bp[,c(1,4)]->mnase_model
  #   colnames(mnase_model)=c("seq","bias")
  #   mnase_vec=mnase_model$bias
  #   names(mnase_vec)=mnase_model$seq
  # }
  data.table::data.table(chr = unique(chromosome), start = (position-(kmer_size/2)+1 ), end = (position+(kmer_size/2))) -> gene_locus_dt_bp
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
