#' Convert InteractionSet to GenomicInteractions Object
#'
#' This function takes an `InteractionSet` object and converts it into a `GenomicInteractions` object. It extracts the first and second anchors from the `InteractionSet` to create the `GenomicInteractions` object.
#'
#' @param iset An `InteractionSet` object to be converted.
#'
#' @return A `GenomicInteractions` object containing the interactions defined in the `iset` object.
#'
#' @examples
#' # Assuming `iset` is an InteractionSet object
#' gint <- iset_to_gint(iset)
#'
#' @importFrom GenomicInteractions GenomicInteractions
#' @importFrom InteractionSet anchors
#' @export
iset_to_gint <- function(iset) {
  GenomicInteractions::GenomicInteractions(anchors(iset)$first, anchors(iset)$second)
}
