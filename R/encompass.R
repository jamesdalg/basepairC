#' Get Encompassing Region of a Genomic Ranges Object
#'
#'@keywords genomicranges iranges
#'@importFrom GenomicRanges seqnames ranges start width GRanges
#'@importFrom IRanges IRanges
#'@param gr A GRanges object
#'@return A single GRange, encompassing the max and min of the input GRanges object
#'@export
encompass=function(gr){
  return(GRanges(seqnames=seqnames(gr)[1],IRanges(gr@ranges@start[1],
                                                  gr@ranges@start[length(gr)]+gr@ranges@width[length(gr)])))
}
