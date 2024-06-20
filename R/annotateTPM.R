annotateTPM <- function(df, org.Mm.eg.db) {
  df$SYMBOL = mapIds(org.Mm.eg.db,
                        keys=row.names(df), 
                        column="SYMBOL",
                        keytype="ENTREZID",
                        multiVals="first")
  SYMBOL <- data.frame(SYMBOL = df$SYMBOL)
  final_df <- cbind(SYMBOL, df[,1:(ncol(df)-1)])
  
  return(final_df)
}