# This script is to generate figures of 16S rRNA sequencing

# Library ====
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)
## Import functions
source(here::here("R", "Functions.R"))

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

# Generate dataframe for plot ====
merged_df_list <- list()
for (sample in sample_list) {
  df <- otu_list[[sample]]
  sampleSheet <- sample_sheet_list[[sample]]
  # Convert cecum_df to a long format
  df_long <- df |>
    tibble::rownames_to_column(var = "otu") |> # Convert row names to a column
    tidyr::pivot_longer(cols = -otu, names_to = "sample", values_to = "abundance")
  
  # Merge with sampleSheet to get age information
  merged_df <- df_long |>
    dplyr::left_join(sampleSheet, by = "sample")
  
  ## Sub the otu description
  if (level == "order") {
    merged_df <- merged_df |>
      mutate(
        otu = sub(".*o__", "", otu),                # Extract text after "f__", "g__" or #o__"
        otu = ifelse(otu == "" | is.na(otu), NA, otu)  # Replace blanks with NA
      )
  } else if (level == "family") {
    merged_df <- merged_df |>
      mutate(
        otu = sub(".*f__", "", otu),                # Extract text after "f__", "g__" or #o__"
        otu = ifelse(otu == "" | is.na(otu), NA, otu)  # Replace blanks with NA
      )
  } else if (level == "genus") {
    merged_df <- merged_df |>
      mutate(
        otu = sub(".*g__", "", otu),                # Extract text after "f__", "g__" or #o__"
        otu = ifelse(otu == "" | is.na(otu), NA, otu)  # Replace blanks with NA
      )
  }
  
  merged_df <- merged_df |>
    group_by(otu) |>
    mutate(
      otu = ifelse(is.na(otu), paste0("Unknown", row_number()), otu)
    ) |>
    ungroup()
  
  ## Calculate total abundance for each OTU across all samples
  top_otu <- merged_df |>
    group_by(otu) |>
    summarize(total_abundance = sum(abundance, na.rm = TRUE)) |>
    arrange(desc(total_abundance)) |>
    slice(1:20) |>
    pull(otu)
  
  ## Filter the original data to include only the top 10 genera
  merged_df <- merged_df |>
    mutate(otu_label = ifelse(otu %in% top_otu, otu, "Other"))
  
  merged_df_list[[sample]] <- merged_df
}



# Plot ====
# Combine OTUs from both datasets to create a consistent color palette
all_otus <- unique(unlist(lapply(merged_df_list, function(df) unique(df$otu_label))))

# Generate a consistent color palette for all OTUs
set.seed(42) # For reproducibility
otu_colors <- setNames(colorRampPalette(brewer.pal(8, "Set3"))(length(all_otus)), all_otus)
for (sample in sample_list) {
  merged_df <- merged_df_list[[sample]]
  
  # Create a stacked bar plot grouped by age
  ggplot(merged_df, aes(x = age, y = abundance, fill = otu_label)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(
      title = paste0("Relative Abundance of ", level, " from ", sample, " by Age"),
      x = "Age",
      y = "Relative Abundance (%)"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10)
    ) +
    guides(fill = guide_legend(override.aes = list(size = 5))) +
    scale_fill_manual(values = otu_colors)  # Apply the consistent color palette
  
  # Save the plot
  ggsave(here::here(figure_dir, paste0("Stack_plot_", sample, "_", level, ".pdf")),
         width = 12, height = 8, dpi = 300)
}