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

## Venn
library(VennDiagram)

## UpSet
library(UpSetR)

## Plot
library(ggplot2)
library(RColorBrewer)
## Import functions
source("/Users/thsieh/TianLab_RNAseq/R/Functions.R")

# Dir ====
doc_dir <- here::here("Doc")
result_dir <- here::here("Results")
figure_dir <- here::here("Results", "Figures")

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

# Filtered count matrix ====
count_filtered <- getFilteredCounts(min_abundance = 5, min_proportion = 0.5, count)

count_annotated <- annotateGENEID(count_filtered, org.Mm.eg.db)

# Perform DEG ====
## Normal aging ====
gf_spf <- "SPF"
for (organ in c("SI", "LI", "Ce")) {
  for (condition in list(c(17, 3), c(78, 3))) {
    sample_sheet_subset <- sampleSheet |>
      dplyr::filter(Organ == !!organ) |>
      dplyr::filter(GF.SPF == !!gf_spf) |>
      dplyr::mutate(condition = Weeks.factor) |>
      dplyr::filter(condition %in% !!condition)
    count_filtered_subset <- count_filtered[,rownames(sample_sheet_subset)]
    dea_results <- performCountDEAByDeseq2(
      sample_sheet = sample_sheet_subset, 
      count = count_filtered_subset, 
      condition = condition,
      org.Mm.eg.db = org.Mm.eg.db,
      result_dir = here::here(
        result_dir, "Statistics", "DEG", 
        paste0(organ, "_", paste0(condition[1], "v", condition[2]), ".csv"))
      )
  }
}

## GF vs SPF ====
condition <- c("GF", "SPF")
for (organ in c("SI", "LI", "Ce")) {
  for (week in c(3, 17, 78)) {
    sample_sheet_subset <- sampleSheet |>
      dplyr::filter(Organ == !!organ) |>
      dplyr::filter(Weeks == !!week) |>
      dplyr::mutate(condition = GF.SPF)
    count_filtered_subset <- count_filtered[,rownames(sample_sheet_subset)]
    dea_results <- performCountDEAByDeseq2(
      sample_sheet = sample_sheet_subset, 
      count = count_filtered_subset, 
      condition = condition,
      org.Mm.eg.db = org.Mm.eg.db,
      result_dir = here::here(
        result_dir, "Statistics", "DEG", 
        paste0(organ, "_", week, "_", paste0(condition[1], "v", condition[2]), ".csv"))
    )
  }
}