# Library ====
library(dplyr)
library(ggplot2)
library(pheatmap)
source(here::here("R", "Functions.R"))

# Dir ====
doc_dir <- here::here("Doc")
gmt_dir <- here::here(doc_dir, "gmt")
result_dir <- here::here("Results")
deg_dir <- here::here(result_dir, "Statistics", "DEG")
fgsea_dir <- here::here(result_dir, "Statistics", "FGSEA")
figure_dir <- here::here(result_dir, "Figures")

# Read FGSEA results ====
suffix <- "DEG_microbiota"
fgsea_list <- list()
for (pathway in c("hallmark", "c2_reactome", "c5_gobp")) {
  for (organ in c("SI", "LI", "Ce")) {
    for (age in c(3, 17, 78)) {
      condition <- paste0(suffix, "_", pathway, "_", organ, "_", age)
      fname <- paste0(condition, "_fgsea_results.csv")
      df <- read.csv(here::here(fgsea_dir, fname))
      fgsea_list[[condition]] <- df
    }
  }
}

# Plot ====
# Heatmap ====
suffix <- "DEG_microbiota"
for (pathway in c("hallmark", "c2_reactome", "c5_gobp")) {
  for (organ in c("SI", "LI", "Ce")) {
    condition <- paste0(suffix, "_", pathway, "_", organ)
    condition_3 <- paste0(condition, "_", 3)
    condition_17 <- paste0(condition, "_", 17)
    condition_78 <- paste0(condition, "_", 78)
    
    # Read the dataframe
    fgsea_3 <- fgsea_list[[condition_3]]
    fgsea_17 <- fgsea_list[[condition_17]]
    fgsea_78 <- fgsea_list[[condition_78]]
    
    # Select top 10 up-regulated pathways (highest NES) & top 10 down-regulated pathways (lowest NES) based on LI_78
    # index pathways
    index <- paste0(suffix, "_", pathway, "_LI_78")
    top_up <- fgsea_list[[index]] %>%
      arrange(padj) %>%
      dplyr::filter(NES > 0) %>%
      slice_head(n = 5)
    
    top_down <- fgsea_list[[index]] %>%
      arrange(padj) %>%
      dplyr::filter(NES < 0) %>%
      slice_head(n = 5)
    
    # Combine selected pathways
    selected_pathways <- c(top_up$pathway, top_down$pathway)
    

    # Filter both dataframes to only include selected pathways
    fgsea_3_filtered <- fgsea_3 %>% filter(pathway %in% selected_pathways)
    fgsea_17_filtered <- fgsea_17 %>% filter(pathway %in% selected_pathways)
    fgsea_78_filtered <- fgsea_78 %>% filter(pathway %in% selected_pathways)
    
    
    # Merge the NES values into a matrix
    heatmap_data <- data.frame(
      pathway = selected_pathways,
      `3` = fgsea_3_filtered$NES[match(selected_pathways, fgsea_3_filtered$pathway)],
      `17` = fgsea_17_filtered$NES[match(selected_pathways, fgsea_17_filtered$pathway)],
      `78` = fgsea_78_filtered$NES[match(selected_pathways, fgsea_78_filtered$pathway)]
    )
    
    colnames(heatmap_data) <- c("pathway", 3, 17, 78)
    
    # Modify row names: remove "GOBP_", "REACTOME_", or "HALLMARK_" and replace "_" with spaces
    cleaned_names <- gsub("^(GOBP_|REACTOME_|HALLMARK_)", "", heatmap_data$pathway)  # Remove prefixes
    cleaned_names <- gsub("_", " ", cleaned_names)  # Replace underscores with spaces
    
    # Assign cleaned names back to rownames
    rownames(heatmap_data) <- cleaned_names
    
    # Define a custom color palette with zero-centered scale
    custom_palette <- colorRampPalette(c("blue", "white", "red"))(50)
    
    heatmap_data <- heatmap_data[,2:4]
    # Generate annotation matrix for p-values
    padj_matrix <- matrix("", nrow = nrow(heatmap_data), ncol = ncol(heatmap_data))
    colnames(padj_matrix) <- colnames(heatmap_data)
    rownames(padj_matrix) <- rownames(heatmap_data)
    
    # Assign significance markers based on padj
    for (i in seq_len(nrow(padj_matrix))) {
      pathway_name <- selected_pathways[i]  # Get pathway name once per row
      
      for (j in seq_len(ncol(padj_matrix))) {
        age_group <- colnames(heatmap_data)[j]  # Get the corresponding age
        
        # Select the correct fgsea result based on age
        fgsea_filtered <- switch(
          age_group,
          "3" = fgsea_3_filtered,
          "17" = fgsea_17_filtered,
          "78" = fgsea_78_filtered
        )
        
        # Extract the padj value
        padj_value <- fgsea_filtered$padj[match(pathway_name, fgsea_filtered$pathway)]
        
        # Assign significance markers based on padj threshold
        if (!is.na(padj_value)) {
          padj_matrix[i, j] <- case_when(
            padj_value < 0.01 ~ "***",
            padj_value < 0.05 ~ "*",
            padj_value < 0.1  ~ "##",
            padj_value < 0.25 ~ "#",
            TRUE              ~ ""  # Default empty string
          )
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
