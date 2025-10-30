freeze=function(melted_df,value_col="value",loc1_col="Var1",loc2_col="Var2"){
  #browser()
    as.formula(paste0(loc1_col,"~",loc2_col)) ->freeze_formula
  reshape2::dcast(melted_df,freeze_formula,value.var=value_col) ->recast_df
  rownames(recast_df) <- recast_df[,loc1_col]

  recast_df[,loc1_col] <- NULL
  recast_df %>% as.matrix() ->recast_matrix
  return(recast_matrix)
}
