#' Get Conditions for Replicate Files
#'
#' This function assigns conditions ("ctrl" for control and "trt" for treatment) to replicate files based on a specified control tag. It returns a list containing the conditions and the ordered replicate file names.
#' The condtiions are ordered by the condition.
#'
#' @param replicate_files A character vector of file names or paths representing the replicate files.
#' @param control_tag A character string representing the tag for control files. Defaults to `"WT"`.
#'
#' @return A list with two elements:
#' \item{conditions}{A character vector representing the conditions for each file ("ctrl" or "trt").}
#' \item{replicate_files}{A character vector of replicate file names, ordered by condition.}
#'
#' @examples
#' replicate_files <- c("WT_sample1.bam", "TRT_sample2.bam", "WT_sample3.bam")
#' result <- get_conditions(replicate_files, control_tag = "WT")
#' print(result)
#'
#' @importFrom data.table data.table
#' @importFrom dplyr arrange
#' @export
get_conditions=function(replicate_files,control_tag="WT"){
  conditions=ifelse(grepl("WT",replicate_files),"ctrl","trt")
  sample_dt=data.table::data.table(file=replicate_files,condition=conditions) %>% dplyr::arrange(condition)
  sample_dt$condition ->conditions
  sample_dt$file->replicate_files
  return(list(conditions=conditions,replicate_files=replicate_files))
}