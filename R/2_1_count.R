# This script is to generate figures of RNA sequencing

# Library ====
library(dplyr)
library(tidyr)
library(tibble)
library(DESeq2)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

## packages for import salmon data
library(tximport)

## import a package for querying SQLite-based annotation data packages
library(AnnotationDbi)

## Link Entrez GeneID to Gene Symbol
library(org.Mm.eg.db)

## Import functions
source("/Users/thsieh/TianLab_RNAseq/R/Functions.R")

# Dir ====
doc_dir <- here::here("Doc")
gtf_dir <- here::here("Doc", "gtf")
result_dir <- here::here("Results")
count_dir <- here::here(result_dir, "Statistics", "RNA_count")

# Doc ====
## Tx info ====
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

## Get keys for each dataset, key is the transcript name (TXNAME)
## Note: select methods will be masked if dplyr is run somewhere
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene <- tx2gene[!is.na(tx2gene$GENEID), ]

## sample sheet ====
sampleSheet <- read.csv(here::here(doc_dir, "./sampleSheet/RNAseq.csv"))
## assign file path
sampleSheet$fpath <- paste(here::here("./Results/Preprocessed/RNAseq/", sampleSheet$Code, "salmon_quant", "quant.sf"))

## sample name
sampleSheet$sample <- paste(sampleSheet$Organ, sampleSheet$GF.SPF, sampleSheet$Weeks, sampleSheet$Replicate, sep = "_")

## new variable
sampleSheet <- sampleSheet |>
  dplyr::mutate(GF.SPF_Weeks = paste0(GF.SPF, "_", Weeks),
                Organ_GF.SPF = paste0(Organ, "_", GF.SPF),
                Weeks.factor = factor(Weeks))

rownames(sampleSheet) <- sampleSheet$sample

## GTF data ====
gene_info <- getGeneInfo(here::here(gtf_dir, "gencode.vM36.basic.annotation.gtf"))

# Count matrix ====
## Check files exist
file.exists(sampleSheet$fpath)

## This package summarize the transcripts into genes
txi <- tximport(sampleSheet$fpath, type = "salmon", tx2gene = tx2gene)

## DESeq2 object
dds <- DESeqDataSetFromTximport(txi,
                                colData = sampleSheet,
                                design = ~ sample)


## Read the count and annotate it
count <- counts(dds)
count <- as.data.frame(count)

colnames(count) <- sampleSheet$sample

count_annotated <- annotateGENEID(count, org.Mm.eg.db)

## Save count
write.csv(count_annotated, here::here(count_dir, "count_all_organ_all_age.csv"), row.names = FALSE)

# Filtered count matrix ====
count_filtered <- getFilteredCounts(min_abundance = 5, min_proportion = 0.5, count)

count_annotated <- annotateGENEID(count_filtered, org.Mm.eg.db)

## Save count
write.csv(count_annotated, here::here(count_dir, "count_filtered_all_organ_all_age.csv"), row.names = FALSE)

# Normalized count matrix ====
dds_normalized <- estimateSizeFactors(dds)  # Ensure size factors are estimated
normalized_count <- counts(dds_normalized, normalized=TRUE)
normalized_count <- as.data.frame(normalized_count)

colnames(normalized_count) <- sampleSheet$sample

normalized_count_annotated <- annotateGENEID(normalized_count, org.Mm.eg.db)

## Save normalized count
write.csv(normalized_count_annotated, here::here(count_dir, "normalized_count_all_organ_all_age.csv"), row.names = FALSE)

# Filtered count matrix ====
count_filtered <- getFilteredCounts(min_abundance = 5, min_proportion = 0.5, normalized_count)

count_annotated <- annotateGENEID(count_filtered, org.Mm.eg.db)

## Save count
write.csv(count_annotated, here::here(count_dir, "normalized_count_filtered_all_organ_all_age.csv"), row.names = FALSE)

# vst-normalized count ====
## Perform variance stabilizing transformation
vst_counts <- vst(dds, blind = TRUE)

## Extract the VST-transformed count matrix
count <- assay(vst_counts)
count <- as.data.frame(count)

## Rename columns using sample names
colnames(count) <- sampleSheet$sample

## Annotate gene IDs
count_annotated <- annotateGENEID(count, org.Mm.eg.db)

## Save VST-normalized count matrix
write.csv(count_annotated, here::here(count_dir, "vst_count_all_organ_all_age.csv"), row.names = FALSE)

# Filtered vst-normalized count matrix ====
count_filtered <- getFilteredCounts(min_abundance = 5, min_proportion = 0.5, count)

count_annotated <- annotateGENEID(count_filtered, org.Mm.eg.db)

## Save count
write.csv(count_annotated, here::here(count_dir, "vst_count_filtered_all_organ_all_age.csv"), row.names = FALSE)