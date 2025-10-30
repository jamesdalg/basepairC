#' Create Data Table with Mappability and GC Content
#'
#' This function creates a data table containing genomic regions specified by the user along with their mappability and GC content. It reads a mappability file for a given genome and calculates mappability and GC content for specified window sizes.
#'
#' @param chr_val The chromosome for which to create the data table (default: "chr15").
#' @param start_val The start position of the genomic region (default: 61982445).
#' @param end_val The end position of the genomic region (default: 61988445).
#' @param genome_val The genome object (e.g., BSgenome object for the specified organism).
#' @param map_fn The file path to the mappability data (default: path to mm10 k24 umap bedgraph).
#' @param mnase_model The file path to the MNase bias model for 6-mers (default: path to mnase_mm10_kmer_map_6mer.txt in basepairC package).
#' @param mnase_model_8bp The file path to the MNase bias model for 8-mers (default: path to mnase_mm10_kmer_map_8mer.txt in basepairC package).
#' @param bsgenome_obj The BSgenome object for the specified organism (default: BSgenome.Mmusculus.UCSC.mm10).
#' 
#' @return A data table with columns for chromosome, start position, end position, sequence, and calculated mappability and GC content for different window sizes.
#'
#' @examples
#' dt_with_map_and_gc <- create_dt_with_map_and_gc(chr_val = "chr15", start_val = 61982445, end_val = 61988445, genome_val = BSgenome.Mmusculus.UCSC.mm10)
#'
#' @importFrom data.table fread as.data.table
#' @importFrom parallel detectCores
#' @importFrom Biostrings getSeq
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom BSgenome.Mmusculus.UCSC.mm10 BSgenome.Mmusculus.UCSC.mm10
#' @importFrom S4Vectors mcols
#' @export
create_dt_with_map_and_gc <- function(chr_val = "chr15", start_val = 61982445, end_val = 61988445, genome_val, map_fn = "/code/R/basepairC/data/mm10/k24.umap.bedgraph",mnase_model=system.file(package="basepairC","extdata/bias_models/mnase_mm10_kmer_map_6mer.txt"),mnase_model_8bp=system.file(package="basepairC","extdata/bias_models/mnase_mm10_kmer_map_8mer.txt"),bsgenome_obj=BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10) {
  print("creating gene locus granges object")
  #browser()
  gene_locus_gr <- data.frame(chr = chr_val, start = start_val, end = end_val) %>% GRanges()
  print("creating gene locus data table")
  gene_locus_gr %>%
    as.data.frame() %>%
    data.table::as.data.table() -> gene_locus_dt
  print("reading mappability data table")
  mappability_dt <- data.table::fread(map_fn, col.names = c("chr", "start", "end", "map"), header = F, nThread = parallel::detectCores())
  print("subsetting mappability data table")
  data.table::setDTthreads(parallel::detectCores())
  mappability_dt[(mappability_dt$chr == chr_val) & (mappability_dt$start >= (start_val - 1e8)) & (mappability_dt$end <= (end_val + 1e8)), ] -> mappability_dt
  print("constructing mappability GRanges object...")
  mappability_gr <- mappability_dt %>% GRanges()
  print("subsetting mappability to gene region")
  mappability_gr %>% subsetByOverlaps(., gene_locus_gr) -> gene_map_gr
  Biostrings::getSeq(bsgenome_obj, gene_locus_gr) -> gene_locus_seq
  data.table::data.table(chr = gene_locus_dt$seqnames %>% unique(), start = seq(gene_locus_dt$start, gene_locus_dt$end), end = seq(gene_locus_dt$start, gene_locus_dt$end), seq = gene_locus_seq %>% as.character() %>% strsplit("") %>% unlist()) -> gene_locus_dt_bp
  
  gene_locus_dt_bp$map_25bp <- map_win(gene_locus_dt_bp$chr, gene_locus_dt_bp$start, window_size = 25, map_gr = mappability_gr)
  gene_locus_dt_bp$map_1bp <- map_win(gene_locus_dt_bp$chr, gene_locus_dt_bp$start, window_size = 1, map_gr = mappability_gr)
  gene_locus_dt_bp$map_1000bp <- map_win(gene_locus_dt_bp$chr, gene_locus_dt_bp$start, window_size = 1000, map_gr = mappability_gr)
  gene_locus_dt_bp$gc_25bp <- gc_win(gene_locus_dt_bp$chr, gene_locus_dt_bp$start, window_size = 25,genome=bsgenome_obj)
  #browser()
  gene_locus_dt_bp$mnase_bias_6bp_viestra <- mnase_bias_6bp(gene_locus_dt_bp$chr, gene_locus_dt_bp$start, window_size = 6,genome=bsgenome_obj) #this is for a model that is on the third nucleotide (2 before, 3 after).
  gene_locus_dt_bp$mnase_bias_6bp <- mnase_bias_kbp(gene_locus_dt_bp$chr, gene_locus_dt_bp$start, kmer_size = 6,mnase_model = mnase_model,genome=bsgenome_obj) #this is for a model that is on the third nucleotide (2 before, 3 after).
  gene_locus_dt_bp$mnase_bias_8bp <- mnase_bias_kbp(gene_locus_dt_bp$chr, gene_locus_dt_bp$start, kmer_size = 8,mnase_model = mnase_model_8bp,genome=bsgenome_obj) #this is for a model that is on the third nucleotide (2 before, 3 after).
  #gene_locus_dt_bp$mnase_bias_6bp <- mnase_bias_6bp(gene_locus_dt_bp$chr, gene_locus_dt_bp$start, window_size = 6)
  #gene_locus_dt_bp$mnase_bias_kmer_6_3_3 <- mnase_bias_kbp(gene_locus_dt_bp$chr, gene_locus_dt_bp$start, kmer_size = 6,mnase_model=mnase_model)
  #gene_locus_dt_bp$gc_1bp=ifelse(dt_with_map_and_gc$seq %in% c("G","C"),1,0)
  #browser()
  gene_locus_dt_bp$gc_1bp=ifelse(gene_locus_dt_bp$seq %in% c("G","C"),1,0)
  return(gene_locus_dt_bp)
}
