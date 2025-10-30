#' Mirror Output Summary
#'
#' This function takes a data frame with genomic coordinates (start1, end1, start2, end2) and creates a mirrored version by swapping the start and end positions between the two sets of coordinates. It then combines the original and mirrored data frames.
#'
#' @param output_summary A data frame with columns `start1`, `end1`, `start2`, and `end2` representing genomic coordinates.
#'
#' @return A data frame that combines the original and mirrored coordinates. The returned data frame is twice the size of the input, with the second half containing the mirrored coordinates.
#'
#' @examples
#' # Assuming df is a data frame with columns start1, end1, start2, end2
#' mirrored_df <- mirror_output_summary(df)
#'
#' @importFrom dplyr mutate select
#' @export
mirror_output_summary <- function(output_summary) {
  rbind(output_summary, output_summary %>% dplyr::mutate(new_start1 = start2, new_end1 = end2, new_start2 = start1, new_end2 = end1) %>% dplyr::mutate(start1 = new_start1, end1 = new_end1, start2 = new_start2, end2 = new_end2) %>% dplyr::select(-new_start1, -new_end1, -new_start2, -new_end2))
}
