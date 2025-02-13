# This script is to provide an overlapping analysis for the significantly enriched GO BP term among
# microbiota-associated DEGs (GF vs SPF, 78w) and aging-associated DEGs (significantly changed slopes 3-17-78w, SPF)

# Library ====
library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gridExtra)
source(here::here("R", "Functions.R"))

# Dir ====
doc_dir <- here::here("Doc")
gmt_dir <- here::here(doc_dir, "gmt")
result_dir <- here::here("Results")
deg_dir <- here::here(result_dir, "Statistics", "DEG")
fgsea_dir <- here::here(result_dir, "Statistics", "FGSEA")
figure_dir <- here::here(result_dir, "Figures")

# Read FGSEA results ====
## FGSEA based on microbiota-associated DEGs (GF vs SPF, 78w) ====
suffix <- "DEG_microbiota"
fgsea_list <- list()
for (pathway in c("hallmark", "c2_reactome", "c5_gobp")) {
  for (organ in c("SI", "LI", "Ce")) {
    for (age in c(78)) {
      condition <- paste0(suffix, "_", pathway, "_", organ, "_", age)
      fname <- paste0(condition, "_fgsea_results.csv")
      df <- read.csv(here::here(fgsea_dir, fname))
      fgsea_list[[condition]] <- df
    }
  }
}

mfgsea_list <- fgsea_list

## FGSEA based on aging trend DEG (3-17-78, SPF) ====
suffix <- "DEG_normal_aging"
fgsea_list <- list()
for (pathway in c("hallmark", "c2_reactome", "c5_gobp")) {
  for (organ in c("SI", "LI", "Ce")) {
    for (age in c("78v3")) {
      condition <- paste0(suffix, "_", pathway, "_", organ, "_", age)
      fname <- paste0(condition, "_fgsea_results.csv")
      df <- read.csv(here::here(fgsea_dir, fname))
      fgsea_list[[condition]] <- df
    }
  }
}

agingfgsea_list <- fgsea_list

# Plot ====

for (pathway in c("hallmark", "c2_reactome", "c5_gobp")) {
  for (organ in c("SI", "LI", "Ce")) {
    
    # Check the condition name
    condition_name1 <- paste0("DEG_normal_aging_", pathway, "_", organ, "_78v3")
    condition_name2 <- paste0("DEG_microbiota_", pathway, "_", organ, "_78")
    
    # Get the fgsea sets
    fgsea_df1 <- agingfgsea_list[[condition_name1]]
    fgsea_df2 <- mfgsea_list[[condition_name2]]
    
    # Set the threshold
    padj_threshold <- 0.05
    
    # Filter significant pathways in both datasets
    significant_fgsea_1 <- fgsea_df1 %>%
      filter(padj < padj_threshold)
    
    significant_fgsea_2 <- fgsea_df2 %>%
      filter(padj < padj_threshold)
    
    # Identify pathways with opposite directions
    significant_fgsea_1_down <- significant_fgsea_1 |>
      filter(NES < 0) |>
      pull(pathway)
    
    significant_fgsea_1_up <- significant_fgsea_1 |>
      filter(NES > 0) |>
      pull(pathway)
    
    significant_fgsea_2_down <- significant_fgsea_2 |>
      filter(NES < 0) |>
      pull(pathway)
    
    significant_fgsea_2_up <- significant_fgsea_2 |>
      filter(NES > 0) |>
      pull(pathway)
    
    # Create named lists for Venn diagram input
    venn_list_1 <- list("Down in d1" = significant_fgsea_1_down, "Up in d2" = significant_fgsea_2_up)
    venn_list_2 <- list("Up in d1" = significant_fgsea_1_up, "Down in d2" = significant_fgsea_2_down)
    
    # Generate the first Venn diagram (Down in D1, Up in D2)
    venn_plot_1 <- ggVennDiagram(venn_list_1, label = "count") +
      ggtitle("Down in aging, Up in GFvsSPF") +
      scale_fill_gradient(low = "white", high = "blue", limits = c(0,60)) +
      theme_classic()
    
    # Generate the second Venn diagram (Up in D1, Down in D2)
    venn_plot_2 <- ggVennDiagram(venn_list_2, label = "count") +
      ggtitle("Up in aging, Down in GFvsSPF") +
      scale_fill_gradient(low = "white", high = "red", limits = c(0,250)) +
      theme_classic()
    
    # Arrange the two plots side by side
    grid.arrange(venn_plot_1, venn_plot_2, ncol = 2)
    
    # Define output file name
    pdf_filename <- here::here(figure_dir, paste0("VennPlot_aging_microbiota_fgsea_", pathway, "_", organ, ".pdf"))
    
    # Save as PDF
    pdf(pdf_filename, width = 12, height = 6)  # Adjust size as needed
    grid.arrange(venn_plot_1, venn_plot_2, ncol = 2)
    dev.off()  # Close the PDF device
  }
}
  
