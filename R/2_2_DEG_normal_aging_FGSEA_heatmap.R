# Library ====
library(dplyr)
library(ggplot2)
library(pheatmap)

# Dir ====
doc_dir <- here::here("Doc")
result_dir <- here::here("Results")
deg_normal_aging_dir <- here::here(result_dir, "Statistics", "DEG_normal_aging")
fgsea_dir <- here::here(result_dir, "Statistics", "FGSEA")
figure_dir <- here::here("Results", "Figures")

# Read FGSEA results ====
suffix <- "DEG_normal_aging"
fgsea_list <- list()
for (pathway in c("hallmark", "c2_reactome", "c5_gobp")) {
  for (organ in c("SI", "LI", "Ce")) {
    for (age in c(17, 78)) {
      condition <- paste0(suffix, "_", pathway, "_", organ, "_", age, "v3")
      fname <- paste0(condition, "_fgsea_results.csv")
      df <- read.csv(here::here(fgsea_dir, fname))
      fgsea_list[[condition]] <- df
    }
  }
}

# Heatmap ====
suffix <- "DEG_normal_aging"
for (pathway in c("hallmark", "c2_reactome", "c5_gobp")) {
  for (organ in c("SI", "LI", "Ce")) {
    condition <- paste0(suffix, "_", pathway, "_", organ)
    condition_17 <- paste0(condition, "_", 17, "v3")
    condition_78 <- paste0(condition, "_", 78, "v3")
    
    # Read the dataframe
    fgsea_17 <- fgsea_list[[condition_17]]
    fgsea_78 <- fgsea_list[[condition_78]]
    
    # Select top 10 up-regulated pathways (highest NES) & top 10 down-regulated pathways (lowest NES)
    top_up <- fgsea_78 %>%
      arrange(desc(NES)) %>%
      slice_head(n = 5)
    
    top_down <- fgsea_78 %>%
      arrange(NES) %>%
      slice_head(n = 5)
    
    # Combine selected pathways
    selected_pathways <- c(top_up$pathway, top_down$pathway)
    
    
    # Filter both dataframes to only include selected pathways
    fgsea_78_filtered <- fgsea_78 %>% filter(pathway %in% selected_pathways)
    fgsea_17_filtered <- fgsea_17 %>% filter(pathway %in% selected_pathways)
    
    # Merge the NES values into a matrix
    heatmap_data <- data.frame(
      Pathway = selected_pathways,
      `17` = fgsea_17_filtered$NES[match(selected_pathways, fgsea_17_filtered$pathway)],
      `78` = fgsea_78_filtered$NES[match(selected_pathways, fgsea_78_filtered$pathway)]
    )
    
    # Modify row names: remove "GO_BP_" and replace "_" with spaces
    cleaned_names <- gsub("^GOBP_", "", heatmap_data$Pathway)  # Remove prefix
    cleaned_names <- gsub("_", " ", cleaned_names)  # Replace underscores with spaces
    
    # Assign cleaned names back to rownames
    rownames(heatmap_data) <- cleaned_names
    
    # Define a custom color palette with zero-centered scale
    custom_palette <- colorRampPalette(c("blue", "white", "red"))(50)
    
    heatmap_data <- heatmap_data[,2:3]
    # Generate annotation matrix for p-values
    padj_matrix <- matrix("", nrow = nrow(heatmap_data), ncol = ncol(heatmap_data))
    colnames(padj_matrix) <- colnames(heatmap_data)
    rownames(padj_matrix) <- rownames(heatmap_data)
    
    # Assign significance markers based on padj
    for (i in 1:nrow(padj_matrix)) {
      for (j in 1:ncol(padj_matrix)) {
        pathway_name <- selected_pathways[i]
        age_group <- colnames(heatmap_data)[j]
        
        # Extract padj value
        padj_value <- ifelse(
          age_group == "17",
          fgsea_17_filtered$padj[fgsea_17_filtered$pathway == pathway_name],
          fgsea_78_filtered$padj[fgsea_78_filtered$pathway == pathway_name]
        )
        
        # Assign markers based on padj threshold
        if (!is.na(padj_value)) {
          if (padj_value < 0.01) {
            padj_matrix[i, j] <- "***"
          } else if (padj_value < 0.05) {
            padj_matrix[i, j] <- "*"
          } else if (padj_value < 0.1) {
            padj_matrix[i, j] <- "##"
          } else if (padj_value < 0.25) {
            padj_matrix[i, j] <- "#"
          }
        }
      }
    }
    
    # Define a custom color palette with zero-centered scale
    custom_palette <- colorRampPalette(c("blue", "white", "red"))(50)
    
    # Ensure consistent breaks centered at 0
    breaks_seq <- seq(min(heatmap_data, na.rm = TRUE), max(heatmap_data, na.rm = TRUE), length.out = 51)
    
    # Define PDF filename
    pdf_filename <- paste0("Heatmap_FGSEA_", suffix, "_", pathway, "_", organ, ".pdf")
    pdf_path <- file.path(figure_dir, pdf_filename)
    
    # Save the heatmap as a PDF
    pdf(pdf_path, width = 12, height = 16)
    pheatmap(
      heatmap_data,
      cluster_rows = FALSE,  # Keep pathway order
      cluster_cols = FALSE,  # Keep sample order fixed
      main = paste0("Top Pathways (", pathway, " - ", organ, ")"),
      color = custom_palette,  # Ensures blue for negative NES and red for positive NES
      breaks = breaks_seq,  # Keeps color mapping consistent
      display_numbers = padj_matrix,  # Show significance labels
      fontsize_number = 10  # Adjust label size
    )
    dev.off()  # Close the PDF device
  }
}
