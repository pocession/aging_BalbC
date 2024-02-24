# This script is to demonstrate how to use GenomicRanges packages to annotate sequence data
# using GenomicRanges

# Library ====
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Transcript data frames ====
transcript_df <- as.data.frame(transcripts(TxDb))
gene_df <- as.data.frame(genes(TxDb))

# Dir ====
inputDir <- here::here("Results/Statistics/DMR")
DMR_list <- list.files(inputDir)
DMR_list <- paste0(inputDir, "/", DMR_list)

# DMR file ====
DMR_file <- DMR_list[[1]]
df <- read.csv(DMR_file)
df$dmr_id <- row.names(df)

# Create GRange objects ====
gr <- GRanges(
  seqnames = df$chr,
  ranges = IRanges(start = df$start, end = df$end),
  pvalue = df$pvalue,
  qvalue = df$qvalue,
  meth.diff = df$meth.diff,
  dmr_id = df$dmr_id
)

# Finding overlap ====
features <- transcripts(TxDb)

# Select the id of closet genes ====
## nearest from peak start
ps.idx <- follow(gr, features)

## nearest from peak end
pe.idx <- precede(gr, features)

## check if any peak not matched to the nearest gene
na.idx <- is.na(ps.idx) & is.na(pe.idx)

## Remove na 
if (sum(na.idx) > 0) { ## suggested by Thomas Schwarzl
  ps.idx <- ps.idx[!na.idx]
  pe.idx <- pe.idx[!na.idx]
  gr <- gr[!na.idx]
}


# Generate the final data frame ====
## Identify the closet gene information
features_df <- data.frame(txStart = start(features[ps.idx]),
                          txEnd = end(features[pe.idx]),
                          txStrand = strand(features[ps.idx]),
                          tx_id = features[ps.idx]$tx_id,
                          tx_name = features[ps.idx]$tx_name)

features_df <- features_df |> 
  dplyr::mutate(distanceToTSS = ifelse(txStrand == "+", txEnd - txStart, txStart - txEnd))

gr_df <- as.data.frame(gr)

## Add the original strand information back
gr_df_tmp <- cbind(gr_df[,1:3], strand = df$strand)
gr_df_tmp <- cbind(gr_df_tmp, gr_df[,6:9])
gr_df <- gr_df_tmp
remove(gr_df_tmp)

final_df <- cbind(gr_df, features_df)
