#' Quick Heatmap Generation
#'
#' This function generates a heatmap from a given matrix using the ComplexHeatmap package. 
#' If the input is not a matrix, it will be converted to a matrix before generating the heatmap.
#' The heatmap does not include clustering or dendrograms, and does not display row or column names.
#'
#' @param mat A matrix or data that can be converted into a matrix.
#' @param direction up or down
#' @param title A title for the heatmap.
#' @param legend_name A name for the legend.
#' @return A heatmap of the input data.
#'
#' @examples
#' mat <- matrix(rnorm(100), 10, 10)
#' quickheat(mat)
#'
#' @importFrom ComplexHeatmap Heatmap
#' @export
quickheat <- function(mat,direction="up",title="",legend_name="value"){
  if(direction=="up"){
    mat=mat[nrow(mat):1,]
  }
  if(any(class(mat) %in% "matrix")) {
    ComplexHeatmap::Heatmap(mat, cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend = F, show_row_names = F, show_column_names = F,name=legend_name,column_title = title,heatmap_legend_param = list(direction = "horizontal") )->hmap
  } else{
    ComplexHeatmap::Heatmap(as.matrix(mat), cluster_rows = F, cluster_columns = F, show_row_dend = F, show_column_dend = F, show_row_names = F, show_column_names = F,name=legend_name,column_title = title,heatmap_legend_param = list(direction = "horizontal") )->hmap
  }
  ComplexHeatmap::draw(hmap, heatmap_legend_side = "bottom")
}
#' Quick Heatmap Generation (ggplot)
#'
#' This function generates a heatmap from a given matrix using the ComplexHeatmap package. 
#' If the input is not a matrix, it will be converted to a matrix before generating the heatmap.
#' The heatmap does not include clustering or dendrograms, and does not display row or column names.
#'
#' @param mat A matrix or data that can be converted into a matrix.
#' @param direction up or down
#' @return A heatmap of the input data.
#'
#' @examples
#' mat <- matrix(rnorm(100), 10, 10)
#' quickheat(mat)
#'
#' @importFrom ComplexHeatmap Heatmap
#' @export
quickheat_gg=function(mat){
  mat %>% reshape2::melt() %>% dplyr::rename(loc1=Var1,loc2=Var2) %>% tidyr::separate(loc1,sep="_",into=c("chr1","start1","end1")) %>% tidyr::separate(loc2,sep="_",into=c("chr2","start2","end2")) %>% dplyr::mutate(start1=as.numeric(start1),start2=as.numeric(start2),end1=as.numeric(end1),end2=as.numeric(end2)) ->mat_long
  plot(  mat_long %>% ggplot(aes(x=start1,y=start2,fill=value)) + geom_raster() + scale_fill_viridis_c( option = "magma", name = "rank(mnase)",direction=1) + theme_minimal() + coord_fixed() + theme(legend.position = "bottom",legend.text = element_text(angle = 90, vjust = 0.5, hjust = 1)) + labs(title="MNase bias matrix, smoothed using a 10bp window",x="Start of bin 1",y="Start of bin 2")
  )
}