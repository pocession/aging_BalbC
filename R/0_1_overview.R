# This script is to provide an overview statistics of the dataset

# Library ====
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

## packages for import salmon data
library(tximport)

## Link Entrez GeneID to Gene Symbol
library(org.Mm.eg.db)

## import packages for differential expression analysis
library(DESeq2)
# Dir ====
doc_dir <- here::here("Doc")
result_dir <- here::here("Results")

## Import functions
source(here::here("R", "Functions.R"))

# Microbiome ====
cecum <- read.csv(here::here(result_dir, "Preprocessed", "Microbiome", "cecum_otu_genus.csv"))
feces <- read.csv(here::here(result_dir, "Preprocessed", "Microbiome", "feces_otu_genus.csv"))

sampleSheet <- data.frame(sample = colnames(cecum)[2:length(colnames(cecum))])
sampleSheet <- sampleSheet |>
  tidyr::separate(sample, c("id", "age", "replicate"), remove = FALSE)

df <- cecum
rownames(df) <- df[,1]
df <- df[,2:ncol(df)]
filtered_df <- getFilteredCounts(min_abundance = 0.01, min_proportion = 0.5, df)
filtered_cecum <- filtered_df

df <- feces
rownames(df) <- df[,1]
df <- df[,2:ncol(df)]
filtered_df <- getFilteredCounts(min_abundance = 0.01, min_proportion = 0.5, df)
filtered_feces <- filtered_df

total_genus <- union(rownames(filtered_cecum), rownames(filtered_feces))

# RNAseq ====
## Annotation objects ====
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

## Get keys for each dataset, key is the transcript name (TXNAME)
## Note: select methods will be masked if dplyr is run somewhere
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2gene <- tx2gene[!is.na(tx2gene$GENEID), ]

## Sample sheet ====
sampleSheet <- read.csv(here::here(doc_dir, "sampleSheet/RNAseq.csv"))

# Arrange the sample sheet
# assign file path
sampleSheet$fpath <- paste(here::here(result_dir, "Preprocessed/RNAseq/", sampleSheet$Code, "salmon_quant", "quant.sf"))

# new variable
sampleSheet <- sampleSheet |>
  dplyr::mutate(GF.SPF_Weeks = paste0(GF.SPF, "_", Weeks),
                Organ_GF.SPF = paste0(Organ, "_", GF.SPF),
                Weeks.factor = factor(Weeks))

# sample name
sampleSheet$sample <- paste(sampleSheet$Organ, sampleSheet$GF.SPF, sampleSheet$Weeks, sampleSheet$Replicate, sep = "_")

rownames(sampleSheet) <- sampleSheet$sample

## Count matrix ====
# Check files exist
file.exists(sampleSheet$fpath)

# This package summarize the transcripts into genes
txi <- tximport(sampleSheet$fpath, type = "salmon", tx2gene = tx2gene)
names(txi$counts) <- sampleSheet$sample

# DESeq2 object
dds <- DESeqDataSetFromTximport(txi,
                                   colData = sampleSheet,
                                   design = ~ sample)


# Read the count and annotate it
dds <- estimateSizeFactors(dds)
count <- counts(dds, normalized=TRUE)
count <- as.data.frame(count)

# Filter out genes with zero counts across all samples
count_filtered <- getFilteredCounts(count, min_abundance = 5, min_proportion = 0.5)

count_annotated <- annotateGENEID(count_filtered, org.Mm.eg.db)

# Save count
write.csv(count, here::here(result_dir, "Statistics", "TPM_all_organ_all_age.csv"), row.names = FALSE)

# RRBS ====
## The methylation count process has been done. 
## It is a long process. We won't do it again.
## Please refer to:
## "Doc/getUnadjustedMethylationCount2018.Rmd" and "Doc/getUnadjustedMethylationCount2022.Rmd"
## "Doc/getAdjustedMethylationCount2018.Rmd" and "Doc/getAdjustedMethylationCount2022.Rmd"

file_list <- list.files(here::here(result_dir, "Statistics", "5hmc_5mc"))
df_list <- list()
for (i in 1:length(file_list)) {
  # Use gsub to remove the suffixes
  cleaned_filenames <- gsub("(_-_1_total.csv|_-_2_total.csv)$", "", file_list[i])
  df <- read.csv(here::here(result_dir, "Statistics", "5hmc_5mc", file_list[i]))
  df <- df |>
    dplyr::mutate(!!cleaned_filenames := freqC) |>
    dplyr::select(chr, start, end, !!cleaned_filenames)

  df_list[[cleaned_filenames]] <- df
}

# Combine all data frames in df_list by 'chr', 'start', and 'end' columns
combined_df <- purrr::reduce(df_list, function(x, y) dplyr::full_join(x, y, by = c("chr", "start", "end")))

# Save methylation count
write.csv(combined_df, here::here(result_dir, "Statistics", "Methylation_all_organ_all_age.csv"), row.names = FALSE)

