# Load necessary libraries ====
library(dplyr)
library(ggplot2)
library(fgsea)
source(here::here("R", "Functions.R"))

# Dir ====
doc_dir <- here::here("Doc")
gmt_dir <- here::here(doc_dir, "gmt")
result_dir <- here::here("Results")
deg_dir <- here::here(result_dir, "Statistics", "DEG")
fgsea_dir <- here::here(result_dir, "Statistics", "FGSEA")
figure_dir <- here::here(result_dir, "Figures")

# Pathways ====
hallmark_pathways <- fgsea::gmtPathways(here::here(gmt_dir, "mh.all.v2024.1.Mm.symbols.gmt"))
c2_reactome_pathways <- fgsea::gmtPathways(here::here(gmt_dir, "m2.cp.reactome.v2024.1.Mm.symbols.gmt"))
c5_gobp_pathways <- fgsea::gmtPathways(here::here(gmt_dir, "m5.go.bp.v2024.1.Mm.symbols.gmt"))
## c3 GTRD (TF)
c3_gtrd_pathways <- fgsea::gmtPathways(here::here(gmt_dir, "m3.gtrd.v2024.1.Mm.symbols.gmt"))

## Define a named list of pathway sets
pathway_sets <- list(
  "hallmark" = hallmark_pathways,
  "c2_reactome" = c2_reactome_pathways,
  "c3_gtrd" = c3_gtrd_pathways,
  "c5_gobp" = c5_gobp_pathways
)

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

# FGSEA ====
## Ranked list ====
## Rank the genes based on -log10(padj)
## Give the positive and negative sign for log2FC > 0 and <0, respectively
rank_list <- list()
deg_list <- deg_list
for (i in 1:length(deg_list)) {
  df <- deg_list[[i]] |>
    # Filter out rows where padj is NA
    dplyr::filter(!is.na(padj)) |>
    # Calculate -log10(padj) and assign positive/negative sign based on log2FoldChange
    dplyr::mutate(rank_value = -log10(padj) * sign(log2FoldChange)) |>
    # Add a small random noise (jitter) to the rank_value to break ties
    dplyr::mutate(rank_value = rank_value + rnorm(n(), sd = 1e-6)) |>
    # Remove non-finite values (NA, NaN, Inf) from the rank_value column
    dplyr::filter(is.finite(rank_value)) |>
    # Arrange the data in descending order of rank_value
    dplyr::arrange(desc(rank_value))
  # Extract the rank_value and assign gene names as names
  rank <- df |>
    dplyr::select(rank_value) |>
    unlist()
  names(rank) <- df$SYMBOL
  rank_list[[names(deg_list)[i]]] <- rank
}

## Run FGSEA ====
suffix <- "DEG_microbiota"
fgsea_results_list <- list()
# Loop through each pathway set
for (pathway_name in names(pathway_sets)) {
  pathways <- pathway_sets[[pathway_name]]
  
  # Loop through each ranked list
  for (i in seq_along(rank_list)) {
    sample_name <- names(rank_list)[i]
    
    # Run fgsea
    fgsea_results <- fgsea(pathways = pathways, 
                           stats = rank_list[[sample_name]],
                           minSize = 15, 
                           maxSize = 500)
    
    # Save results to the list
    fgsea_results_list[[paste0(pathway_name, "_", sample_name)]] <- fgsea_results
    
    # Convert list columns (e.g., leadingEdge) to comma-separated strings
    fgsea_results_cleaned <- fgsea_results
    if ("leadingEdge" %in% colnames(fgsea_results_cleaned)) {
      fgsea_results_cleaned$leadingEdge <- sapply(fgsea_results_cleaned$leadingEdge, paste, collapse = ", ")
    }
    
    # Construct file name with pathway name first
    file_name <- paste0(suffix, "_", pathway_name, "_", sample_name, "_fgsea_results.csv")
    file_path <- file.path(fgsea_dir, file_name)
    
    # Save results as CSV
    write.csv(fgsea_results_cleaned, file = file_path, row.names = FALSE)
  }
}