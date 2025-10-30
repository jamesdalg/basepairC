#' Calculate GC Content in a Sliding Window
#'
#' This function calculates the GC content within a sliding window across specified genomic positions on a chromosome. It is designed to work with genomic data, utilizing a specified genome object for sequence retrieval.
#'
#' @param chromosome The chromosome for which to calculate GC content.
#' @param position A vector of positions on the chromosome to calculate GC content for.
#' @param genome The genome object from which to retrieve sequence data. By default, this is set to the mouse genome (BSgenome.Mmusculus.UCSC.mm10).
#' @param window_size The size of the sliding window for GC content calculation. Default is 5.
#' @param n_cores The number of cores to use for parallel processing. Default is half of the available cores detected on the system.
#'
#' @return A vector of GC content values for each position, calculated within the specified sliding window.
#'
#' @examples
#' # Calculate GC content for positions 100 to 200 on chromosome 1
#' gc_content <- gc_win(
#'   chromosome = "chr1",
#'   position = 100:200,
#'   genome = BSgenome.Mmusculus.UCSC.mm10,
#'   window_size = 5
#' )
#'
#' @importFrom Biostrings getSeq
#' @importFrom data.table data.table
#' @importFrom parallel detectCores
#' @importFrom GenomicRanges GRanges
#' @importFrom pbmcapply pbmclapply
#' @export
gc_win <- function(chromosome, position, genome = BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10, window_size = 5, n_cores = parallel::detectCores() / 2) {
  data.table::data.table(chr = unique(chromosome), start = min(position), end = max(position)) %>% GRanges() -> gene_locus_gr
  Biostrings::getSeq(genome, gene_locus_gr) -> gene_locus_seq
  data.table::data.table(chr = chromosome, start = seq(min(position), max(position))) %>%
    dplyr::mutate(seq = gene_locus_seq %>% as.character() %>% strsplit("") %>% unlist(), end = start) %>%
    dplyr::select(chr, start, end, seq) -> gene_locus_dt_bp
  # browser()
  pbmcapply::pbmclapply(1:nrow(gene_locus_dt_bp), function(i) {
    # print(i)
    min_window_pos <- max(i, i - round(window_size / 2))
    max_window_pos <- min(i + round(window_size / 2), nrow(gene_locus_dt_bp))
    gene_locus_dt_bp[min_window_pos:max_window_pos, "seq"] %>% unlist() -> window_seq
    sum(window_seq %in% c("G", "C")) / length(window_seq) -> gc_content
    return(gc_content)
  }, mc.cores = n_cores) %>% unlist() -> gc_content_window
  return(gc_content_window)
}
