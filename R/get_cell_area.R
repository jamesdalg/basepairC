#' Calculate Area of a Cell Based on ID
#'
#' This function calculates the area of a cell given a cell ID that encodes position data.
#' The cell ID should contain two position data points separated by a space. Each position data
#' is converted to a GRanges object using `CNVScope::underscored_pos_to_GRanges`, from which
#' the width is extracted. The area is computed as the product of the widths of these two
#' positions.
#'
#' @param cell_id A character string representing the cell ID, which includes two position
#'                data points separated by a space. Each position should be in a format
#'                compatible with `CNVScope::underscored_pos_to_GRanges`.
#'
#' @return Returns a numeric value representing the area of the cell.
#'
#' @examples
#' # Example cell ID format: "chr1_100_200 chr2_150_250"
#' get_cell_area("chr1_100_200 chr2_150_250") # Calculates the area based on positions.
#'
#' @importFrom CNVScope underscored_pos_to_GRanges
#' @export
get_cell_area=function(cell_id) {
  (CNVScope::underscored_pos_to_GRanges(strsplit(cell_id," ")[[1]][[1]])@ranges@width)*(CNVScope::underscored_pos_to_GRanges(strsplit(cell_id," ")[[1]][[2]])@ranges@width) ->area
  return(area)
}