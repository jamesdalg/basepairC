#' Normalize Read Counts using edgeR
#'
#' This function takes read count matrices for control and treatment groups,
#' performs normalization using the edgeR package, and returns the TMM-normalized
#' counts for further analysis.
#'
#' @param control A matrix of read counts for the control group.
#' @param treatment A matrix of read counts for the treatment group.
#' @return A list containing two matrices with TMM-normalized read counts for
#'         control and treatment groups.
#' @import edgeR
#' @examples
#' control <- matrix(rpois(200, lambda = 10), ncol = 5)
#' treatment <- matrix(rpois(200, lambda = 20), ncol = 5)
#' result <- normalizeReadCounts(control, treatment)
#' control_normalized <- result$control
#' treatment_normalized <- result$treatment
#' @export
normalizeReadCounts <- function(control, treatment) {
  
  #convert each matrix to a numeric vector
  # Combine matrices and create group factor
  all_counts <- cbind(as.numeric(control), as.numeric(treatment))
  group <- factor(c(rep("Control", ncol(control)), rep("Treatment", ncol(treatment))))

  # Create DGEList object
  dge <- DGEList(counts = all_counts, group = group)

  # Filter lowly expressed genes
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes=FALSE]

  # Normalize for library size using TMM normalization
  dge <- calcNormFactors(dge)

  # Output normalized counts
  normalized_counts <- cpm(dge, log=FALSE)  # Use log=TRUE if you prefer log-transformed outputs

  # Split normalized counts back into control and treatment
  control_normalized <- normalized_counts[, 1:ncol(control)]
  treatment_normalized <- normalized_counts[, (ncol(control) + 1):ncol(all_counts)]
  matrix(control_normalized, nrow = nrow(control), ncol = ncol(control)) -> control_normalized
  matrix(treatment_normalized, nrow = nrow(treatment), ncol = ncol(treatment)) -> treatment_normalized

  # Return as a list of matrices
  return(list(control = control_normalized, treatment = treatment_normalized))
}
