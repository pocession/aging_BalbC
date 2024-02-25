# Annotate DMR regions with mm10

# Library ====
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPseeker)
library(org.Mm.eg.db)

# Assign the mm10 objects ====
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Extract TSS positions ====
tss <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene, columns=c("tx_id", "gene_id"), use.names=TRUE)
tss <- resize(tss, width=1, fix='start')
tss_df <- as.data.frame(tss)
tss_df$transcript <- row.names(tss_df)
colnames(tss_df) <- c("tx_chr","tx_start","tx_end","tx_width","tx_strand","tx_id","gene_id","transcript")

# Dir ====
inputDir <- here::here("Results/Statistics/")

# Read DMR data ====
file_list <- list.files(here::here(inputDir, "DMR"))
file_list <- paste0(inputDir, "/DMR/", file_list)

# Get annotated DMR data frame ====
annotated_DMR_list <- lapply(file_list, .get_annotated_DMR)

.get_annotated_DMR <- function(file) {
  DMR <- read.csv(here::here(file))
  file_name <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file))
  
  ## Filter out non-standard chr
  wanted <- paste0("chr", seq(1:23))
  wanted <- c(wanted, "chrX", "chrY")
  DMR <- DMR |>
    dplyr::filter(chr %in% wanted) |>
    dplyr::mutate(start = as.numeric(start),
                  end = as.numeric(end))
  
  # Save DMR as a bed file
  write.table(DMR, here::here("Results/Tmp", paste0(file_name, "_tmp.bed")), 
              quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  peakfile <- here::here("Results/Tmp", paste0(file_name, "_tmp.bed"))
  peak <- readPeakFile(peakfile)
  peakAnno <- annotatePeak(peakfile,
                           TxDb=txdb, annoDb="org.Mm.eg.db")
  
  annotated_df <- as.data.frame(peakAnno)
  colnames(annotated_df) <- c("dmrChr","dmrStart","dmrEnd","width","strand_none", "dmrStrand","dmrPvalue","dmrQvalue","dmrMethDiff",
                              "annotation","geneChr","geneStart","geneEnd","geneLength","geneStrand","geneId",
                              "transcriptId","distanceToTSS","ENSEMBL","SYMBOL","GENENAME")
  annotated_df <- annotated_df |>
    dplyr::select(-c(strand_none))
  
  
  if (file.exists(peakfile)) {
    #Delete file if it exists
    file.remove(peakfile)
  }
  write.table(annotated_df, here::here("Results/Tmp", paste0(file_name,"_annotated.bed")), 
              quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  return(annotated_df)
}
