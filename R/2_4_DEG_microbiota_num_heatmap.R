# Load necessary libraries ====
library(dplyr)
library(ggplot2)
library(pheatmap)
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

# Create a number of DEGs for each tissue ====
deg_df <- matrix(ncol=5, nrow=0)
colnames(deg_df) <- c("sample", "Tissue", "Weeks", "direction", "num_deg")
for (i in 1:length(deg_filtered_list)) {
  sample <- names(deg_filtered_list[i])
  condition <- unlist(strsplit(names(deg_filtered_list[i]), "_"))
  
  # Get the numbers of DEGs
  num_up_deg <- deg_filtered_list[[i]] |>
    dplyr::filter(log2FoldChange > 0) |>
    nrow()
  
  num_down_deg <- deg_filtered_list[[i]] |>
    dplyr::filter(log2FoldChange < 0) |>
    nrow()
  
  tmp_df <- data.frame(sample =  c(sample, sample),
                       Tissue = c(condition[1], condition[1]),
                       Weeks = c(condition[2], condition[2]),
                       direction =c("up", "down"),
                       num_deg = c(num_up_deg, num_down_deg))
  
  deg_df <- rbind(deg_df, tmp_df)
}
remove(tmp_df, i)

# Annotation of heatmap
annot_col <- data.frame(row.names = sampleSheet$sample,
                        GF.SPF_Weeks = sampleSheet$GF.SPF_Weeks,
                        Organ = sampleSheet$Organ)
# upregulated genes
# Generate matrix and make figures
deg_df_subset <- deg_df |> 
  dplyr::filter(direction == "up") |> 
  dplyr::select(Tissue, Weeks, num_deg)
deg_mtx <- reshape(deg_df_subset, idvar = "Tissue", timevar = "Weeks", direction = "wide")
deg_mtx <- as.matrix(deg_mtx[,2:ncol(deg_mtx)])
colnames(deg_mtx) <- c(3, 17, 78)
rownames(deg_mtx) <- c("SI", "Ce", "LI")

p <- pheatmap(deg_mtx,fontsize = 4,
              colorRampPalette(c("white","red"))(100),
              cluster_cols=FALSE, cluster_rows=FALSE,scale="row")
# cannot directly save the pheatmap object to pdf
savePheatmapPdf(p, here::here(figure_dir, "Heatmap_deg_microbiota_up_num.pdf"))
p

# downregulated genes
# Generate matrix and make figures
deg_df_subset <- deg_df |> 
  dplyr::filter(direction == "down") |> 
  dplyr::select(Tissue, Weeks, num_deg)
deg_mtx <- reshape(deg_df_subset, idvar = "Tissue", timevar = "Weeks", direction = "wide")
deg_mtx <- as.matrix(deg_mtx[,2:ncol(deg_mtx)])
colnames(deg_mtx) <- c(3, 17, 78)
rownames(deg_mtx) <- c("SI", "Ce", "LI")

p <- pheatmap(deg_mtx,fontsize = 4,
              colorRampPalette(c("white","blue"))(100),
              cluster_cols=FALSE, cluster_rows=FALSE,scale="row")
# cannot directly save the pheatmap object to pdf
savePheatmapPdf(p, here::here(figure_dir, "Heatmap_deg_microbiota_down_num.pdf"))
p
