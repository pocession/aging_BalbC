# This script is to generate figures of 16S rRNA sequencing

# Library ====
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
## Import functions
source("~/TianLab_RNAseq/R/Functions.R")

# Dir ====
doc_dir <- here::here("Doc")
result_dir <- here::here("Results")
figure_dir <- here::here("Results", "Figures")

# Doc ====


# Read OTU ====
sample_list <- c("cecum", "feces")
# level <- "order"
level <- "genus"
otu_list <- list()
for (sample in sample_list) {
  otu_df <- read.csv(here::here(result_dir, "Preprocessed", "Microbiome", paste0(sample, "_otu_", level, ".csv")))
  rownames(otu_df) <- otu_df[,1]
  otu_df <- otu_df[,2:ncol(otu_df)]
  otu_list[[sample]] <- otu_df
}

# Generate sample sheet ====
sample_sheet_list <- list()
for (sample in sample_list) {
  df <- otu_list[[sample]]
  sampleSheet <- data.frame(sample = colnames(df))
  sampleSheet <- sampleSheet |>
    tidyr::separate(sample, c("id", "age", "replicate"), remove = FALSE) |>
    dplyr::mutate(Sample.ID = sample)
  sample_sheet_list[[sample]] <- sampleSheet
}


# PCA ====
for (sample in sample_list) {
  df <- otu_list[[sample]]
  sample_sheet <- sample_sheet_list[[sample]]
  df_t <- t(df)
  pca_result <- prcomp(df_t, center = TRUE, scale. = FALSE) # some zero-abundance gene affect the scaling
  
  ## Create a data frame for PCA results
  ## Combine PCA results with sample information for plotting
  pca_data_combined <- getCombinedSamplePCAResults(pca_result, sample_sheet)
  
  ## Plot PCA by Ages
  ## Define the custom shapes and colors for each gender
  custom_shapes <- c("17w" = 15, "18m" = 17, "24m" = 24) #age
  custom_colors <- c("17w" = "blue", "18m" = "black", "24m" = "red") #age
  
  getPCAPlots(
    pca_data_combined = pca_data_combined , 
    pca_result = pca_result,
    figure_dir = here::here(figure_dir),
    output_filename = paste0("PCA_", level, "_", sample, ".pdf"),
    shape_factor = "age",
    color_factor = "age",
    custom_shapes = custom_shapes,
    custom_colors = custom_colors,
    plot_title = paste0("PCA of ", sample, " in ", level, " level by ages"),
    label_samples = FALSE,
    save_plot = TRUE
  )
}