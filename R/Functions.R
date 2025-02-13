#' annotateGENEID 

annotateGENEID <- function(df, org.Mm.eg.db) {
  df$SYMBOL = mapIds(org.Mm.eg.db,
                     keys=row.names(df), 
                     column="SYMBOL",
                     keytype="ENTREZID",
                     multiVals="first")
  SYMBOL <- data.frame(SYMBOL = df$SYMBOL)
  final_df <- cbind(SYMBOL, df[,1:(ncol(df)-1)])
  
  return(final_df)
}

#' getFilteredCounts
#' @description
#' This function generates a filtered count matrix based on the given parameters
#' @parameter min_abundance the minimal abundance to retain a gene, default is 5
#' @parameter min_proportion the percentage of samples where each gene has counts >= min_abundance
#' defaut is 0.5 (50%)
#' @parameter count
#' @return a filtered dataframe. 
#' @example 
#' filtered_count <-  getFilteredCounts(min_abundance = 5, min_proportion = 0.5, count
getFilteredCounts <- function(min_abundance = 5, min_proportion = 0.5, count) {
  
  ## Filter genes based on minimum abundance and proportion
  ## Calculate the number of samples where each gene has counts >= min_abundance
  keep_genes <- rowSums(count >= min_abundance) >= (min_proportion * ncol(count))
  
  ## Subset the count matrix to keep only the genes that meet the criteria
  filtered_count <- count[keep_genes, ]
  return(filtered_count)
}

# savePheatmapPdf
# Define a function to save the pheatmap object to a PDF
savePheatmapPdf <- function(x, filename, width = 7, height = 7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}