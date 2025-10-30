#' Calculate Mappability Content in a Sliding Window
#'
#' This function calculates the mean mappability content within a sliding window for specified genomic positions on a chromosome with a `GenomicRanges` object that specifies mappability.
#'
#' @param chromosome The chromosome for which to calculate mappability content. Defaults to the chromosome column in `gene_locus_dt_bp`.
#' @param position A vector of positions on the chromosome to calculate mappability content for. Defaults to the start column in `gene_locus_dt_bp`.
#' @param map_gr A `GenomicRanges` object containing mappability data for the genome.
#' @param window_size The size of the sliding window for mappability content calculation. Default is 200.
#' @param n_cores The number of cores to use for parallel processing. Default is half of the available cores detected on the system.
#'
#' @return A vector of mappability content values for each position, calculated within the specified sliding window.
#'
#' @examples
#' # Assuming `gene_locus_dt_bp` and `mappability_gr` are defined
#' mappability_content <- map_win(
#'   chromosome = gene_locus_dt_bp$chr,
#'   position = gene_locus_dt_bp$start,
#'   map_gr = mappability_gr,
#'   window_size = 200
#' )
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom data.table data.table
#' @importFrom dplyr mutate select
#' @importFrom pbmcapply pbmclapply
#' @importFrom parallel detectCores
#' @importFrom S4Vectors mcols
#' @export
map_win <- function(chromosome = gene_locus_dt_bp$chr, position = gene_locus_dt_bp$start, map_gr = mappability_gr, window_size = 200, n_cores = parallel::detectCores() / 2) {
  #browser()
  data.table::data.table(chr = unique(chromosome), start = min(position), end = max(position)) %>% GRanges() -> gene_locus_gr
  map_gr %>% subsetByOverlaps(., gene_locus_gr + min(1e5, window_size * 2)) -> map_gr
   

  # Biostrings::getSeq(genome,gene_locus_gr)->gene_locus_seq
  data.table::data.table(chr = chromosome, start = seq(min(position), max(position))) %>%
    dplyr::mutate(end = start) %>%
    dplyr::select(chr, start, end) -> gene_locus_dt_bp

  gene_locus_dt_bp %>% GRanges() -> gene_locus_dt_bp_gr
  # data.table(chr=unique(chromosome),start=min(position),end=max(position)) %>% GRanges()->gene_locus_dt_bp_gr
  # browser()
  pbmcapply::pbmclapply(1:nrow(gene_locus_dt_bp), function(i) {
    # print(i)
    min_window_pos <- max(i, i - round(window_size / 2))
    max_window_pos <- min(i + round(window_size / 2), nrow(gene_locus_dt_bp))
    mean(mcols(subsetByOverlaps(map_gr, gene_locus_dt_bp_gr[min_window_pos:max_window_pos]))$map, na.rm = T) -> mappability
    return(mappability)
  }, mc.cores = n_cores) %>% unlist() -> map_content_window
  return(map_content_window)
}
