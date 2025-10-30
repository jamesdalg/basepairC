#' Convert a Melted Data Frame into a Matrix
#'
#' This function converts a melted data frame into a matrix, where the rows and columns correspond to specified location columns.
#'
#' @param melted_df A data frame in a melted format, typically with columns representing values and their associated locations.
#' @param value_col A character string specifying the column name that contains the values to populate the matrix. Default is `"value"`.
#' @param loc1_col A character string specifying the column name for the row locations of the matrix. Default is `"Var1"`.
#' @param loc2_col A character string specifying the column name for the column locations of the matrix. Default is `"Var2"`.
#'
#' @return A matrix where rows are indexed by the values in `loc1_col`, columns are indexed by the values in `loc2_col`, and the matrix elements are populated by the values in `value_col`.
#'
#' @details 
#' The function takes a melted data frame and casts it into a wide format using the specified location columns. The resulting data frame is then converted into a matrix, with row names set by the values in `loc1_col`.
#'
#' @examples
#' \dontrun{
#' melted_df <- data.frame(Var1 = c("A", "A", "B", "B"), Var2 = c("X", "Y", "X", "Y"), value = 1:4)
#' matrix_result <- freeze(melted_df, value_col = "value", loc1_col = "Var1", loc2_col = "Var2")
#' }
#'
#' @importFrom reshape2 dcast
#' @export
freeze=function(melted_df,value_col="value",loc1_col="Var1",loc2_col="Var2"){
  #browser()
    as.formula(paste0(loc1_col,"~",loc2_col)) ->freeze_formula
  reshape2::dcast(melted_df,freeze_formula,value.var=value_col) ->recast_df
  rownames(recast_df) <- recast_df[,loc1_col]

  recast_df[,loc1_col] <- NULL
  recast_df %>% as.matrix() ->recast_matrix
  return(recast_matrix)
}
