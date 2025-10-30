#' Mix and Combine Segments from Two Objects
#'
#' This function takes two segment objects (lists of vectors), combines the segments element-wise, removes duplicates, and sorts the resulting segments.
#'
#' @param seg_ob1 A list of vectors, where each vector contains segments to be combined.
#' @param seg_ob2 Another list of vectors, of the same length as `seg_ob1`, to be combined with `seg_ob1`.
#'
#' @return A list of vectors, where each vector contains the unique and sorted combination of the corresponding elements from `seg_ob1` and `seg_ob2`.
#'
#' @details
#' The function first determines the minimum length between `seg_ob1` and `seg_ob2` to ensure that elements are combined safely. 
#' For each element pair (from `seg_ob1` and `seg_ob2`), it combines the segments, removes duplicates, and sorts the result in ascending order. 
#' The final result is returned as a list of these combined segments.
#'
#' @importFrom magrittr %>%
#' 
#' @examples
#' seg1 <- list(c(1, 2, 3), c(5, 6))
#' seg2 <- list(c(2, 4), c(6, 7))
#' mixed_segments <- mix_segments(seg1, seg2)
#' # Returns list: list(c(1, 2, 3, 4), c(5, 6, 7))
#'
#' @export
mix_segments=function(seg_ob1,seg_ob2){
  obj_len=min(length(seg_ob1),length(seg_ob2))
  lapply(1:obj_len,function(i){
    seg_combined_obj_item=c(seg_ob1[[i]],seg_ob2[[i]]) %>% unique() %>% sort()
    return(seg_combined_obj_item)
  })->seg_combined_obj
  names(seg_combined_obj)=names(seg_ob1)
  return(seg_combined_obj)
}
# merge_segments=function(seg_ob,min_seg_size=50,matrix_max=6000){
#   obj_len=min(length(seg_ob))
#   lapply(1:obj_len,function(i){
#     diff(seg_ob[[i]])
#     seg_ob[[i]][-which((diff(c(1,seg_ob[[i]],matrix_max)))<min_seg_size)]->seg_ob_merged
#     return(seg_ob_merged)
#   }->seg_ob_merged_combined
#   names(seg_ob_merged_combined)=names(seg_ob)
#   return(seg_ob_merged_combined)
# }
