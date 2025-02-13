# Load necessary libraries ====
library(dplyr)
library(ggplot2)
library(reshape2)
library(UpSetR)
source(here::here("R", "Functions.R"))

# Dir ====
doc_dir <- here::here("Doc")
result_dir <- here::here("Results")
deg_dir <- here::here(result_dir, "Statistics", "DEG")
count_dir <- here::here(result_dir, "Statistics", "RNA_count")
figure_dir <- here::here(result_dir, "Figures")

# sample sheet ====
sampleSheet <- read.csv(here::here(doc_dir, "sampleSheet", "RNAseq.csv"))

# Arrange the sample sheet
# assign file path
sampleSheet$fpath <- paste(here::here("./Results/Preprocessed/RNAseq/", sampleSheet$Code, "salmon_quant", "quant.sf"))

# new variable
sampleSheet <- sampleSheet |>
  dplyr::mutate(GF.SPF_Weeks = paste0(GF.SPF, "_", Weeks),
                Organ_GF.SPF = paste0(Organ, "_", GF.SPF),
                Weeks.factor = factor(Weeks))

# sample name
sampleSheet$sample <- paste(sampleSheet$Organ, sampleSheet$GF.SPF, sampleSheet$Weeks, sampleSheet$Replicate, sep = "_")

# Read microbiota DEG ====
deg_list <- list()
for (tissue in c("SI", "Ce", "LI")) {
  for (age in c(3, 17, 78)) {
    condition <- paste0(tissue, "_", age)
    fname <- paste0(condition, "_GFvSPF.csv")
    df <- read.csv(here::here(deg_dir, fname))
    deg_list[[condition]] <- df
  }
}

# Filter DEG ====
deg_filtered_list <- list()
for (tissue in c("SI", "Ce", "LI")) {
  for (age in c(3, 17, 78)) {
    condition <- paste0(tissue, "_", age)
    df <- deg_list[[condition]]
    filtered_df <- df |>
      dplyr::filter(padj < 0.05) |>
      dplyr::filter(abs(log2FoldChange) > 0.5)
    deg_filtered_list[[condition]] <- filtered_df
  }
}

# Plot ====
gene_sets <- list()
for (i in 1:length(deg_filtered_list)) {
  condition <- names(deg_filtered_list[i])
  gene_sets[[condition]] <- unlist(deg_filtered_list[[i]]$SYMBOL)
}

all_genes <- unique(unlist(gene_sets))
binary_matrix <- sapply(gene_sets, function(x) as.integer(all_genes %in% x))
colnames(binary_matrix) <- names(gene_sets)
rownames(binary_matrix) <- all_genes
binary_matrix <- as.data.frame(binary_matrix)

p <- UpSetR::upset(binary_matrix, sets = names(gene_sets), 
                   order.by = "freq", 
                   main.bar.color = "blue", 
                   sets.bar.color = "red",
                   sets.x.label = "Gene Sets",
                   text.scale = c(2, 1.5, 1, 1, 1.5, 1.2))
p
pdf(here::here(figure_dir, "Upset_deg_microbiota.pdf"), height = 6, width = 6.5)
print(p) # use print here
dev.off()