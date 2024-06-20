annotateDEA <- function(dds, org.Mm.eg.db) {
  res <- results(dds)
  
  res$symbol = mapIds(org.Mm.eg.db,
                      keys=row.names(res), 
                      column="SYMBOL",
                      keytype="ENTREZID",
                      multiVals="first")
  
  res$name =   mapIds(org.Mm.eg.db,
                      keys=row.names(res), 
                      column="GENENAME",
                      keytype="ENTREZID",
                      multiVals="first")
  return(as.data.frame(res))
}