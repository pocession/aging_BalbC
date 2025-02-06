# Load necessary libraries ====
library(dplyr)
library(ggplot2)
library(pheatmap)

# Dir ====
doc_dir <- here::here("Doc")
gmt_dir <- here::here(doc_dir, "gmt")
result_dir <- here::here("Results")
deg_age_trend_dir <- here::here(result_dir, "Statistics", "DEG_age_trend")
count_dir <- here::here(result_dir, "Statistics", "RNA_count")
figure_dir <- here::here(result_dir, "Figures")

# Read the count matrix
file_path <- here::here(count_dir, "normalized_count_filtered_all_organ_all_age.csv")
count_matrix <- read.csv(file_path)

# Remove duplicate gene symbols by keeping the first occurrence
count_matrix <- count_matrix[!duplicated(count_matrix$SYMBOL), ]

rownames(count_matrix) <- count_matrix$SYMBOL
count_matrix <- count_matrix[,2:ncol(count_matrix)]

# Step 1: Subset samples with SPF treatment
spf_samples <- grep("SPF", colnames(count_matrix), value = TRUE)
count_spf <- count_matrix[, spf_samples]

# Step 2: Separate samples into SI, LI, and Ce categories
si_samples <- grep("^SI_", colnames(count_spf), value = TRUE)
li_samples <- grep("^LI_", colnames(count_spf), value = TRUE)
ce_samples <- grep("^Ce_", colnames(count_spf), value = TRUE)

count_si <- count_spf[, si_samples]
count_li <- count_spf[, li_samples]
count_ce <- count_spf[, ce_samples]

# Read DEG data ====
deg_si <- read.csv(here::here(deg_age_trend_dir, "deg_si_age_trend.csv"))
deg_li <- read.csv(here::here(deg_age_trend_dir, "deg_li_age_trend.csv"))
deg_ce <- read.csv(here::here(deg_age_trend_dir, "deg_ce_age_trend.csv"))

# XY scatter ====
# Extract the top up- and down-regulated genes
top_up_gene <- deg_li$Gene[which.max(deg_li$slope)]   # Highest slope
top_down_gene <- deg_li$Gene[which.min(deg_li$slope)] # Lowest slope

# Extract expression data for both genes
top_genes <- c(top_up_gene, top_down_gene)
expression_data <- count_li[rownames(count_li) %in% top_genes, , drop = FALSE]

# Extract ages from column names
ages <- as.numeric(sub(".*_(\\d+)_\\d+$", "\\1", colnames(expression_data)))

# Transform expression data to long format
expression_long <- as.data.frame(t(expression_data))
colnames(expression_long) <- top_genes
expression_long$Age <- ages

# Reshape for ggplot
expression_long <- pivot_longer(expression_long, cols = top_genes, names_to = "Gene", values_to = "Expression")

# Generate the scatter plot with trend lines
p <- ggplot(expression_long, aes(x = Age, y = Expression, color = Gene)) +
  geom_point(size = 3) +  # Scatter points
  geom_smooth(method = "lm", se = FALSE, linetype = "solid", size = 1) +  # Add regression trendline
  scale_y_continuous(trans = "log2") +  # Log2 transform the Y-axis
  labs(
    title = "Top Up- and Down-Regulated Genes in LI Tissue",
    x = "Age",
    y = "log2(Normalized Expression)",
    color = "Gene"
  ) +
  theme_minimal() +
  ylim(0,1024) +
  scale_color_manual(values = c("blue", "red"))  # Blue for down, red for up

# Save the plot as a PDF
pdf(file = here::here(figure_dir, "Scatter_top_deg_age_trend_deg_li.pdf"), width = 6, height = 6)
print(p)
dev.off()

# Piechart ====
# Get unique gene counts across all tissues
unique_genes <- unique(c(deg_si$Gene, deg_li$Gene, deg_ce$Gene))

# Count the number of unique DEGs in each tissue
deg_counts <- c(
  SI = length(intersect(unique_genes, deg_si$Gene)),
  LI = length(intersect(unique_genes, deg_li$Gene)),
  Ce = length(intersect(unique_genes, deg_ce$Gene))
)

# Calculate the total number of unique DEGs
total_unique_genes <- length(unique_genes)

# Convert to dataframe for ggplot
deg_df <- data.frame(
  Tissue = names(deg_counts),
  Count = as.numeric(deg_counts)
)

# Calculate percentages based on unique genes
deg_df$Percentage <- round((deg_df$Count / total_unique_genes) * 100, 1)

# Create label with count and percentage
deg_df$Label <- paste0(deg_df$Tissue, "\n", deg_df$Count, " (", deg_df$Percentage, "%)")

# Generate the pie chart
p <- ggplot(deg_df, aes(x = "", y = Count, fill = Tissue)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar(theta = "y") +
  theme_void() +  # Remove background grid and axis
  labs(title = paste0("Proportion of Unique DEGs Across Tissues\nTotal: ", total_unique_genes, " Genes")) + 
  geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 5) +
  scale_fill_brewer(palette = "Set2")

# Save the plot as a PDF
pdf(file = here::here(figure_dir, "Piechart_normal_aging_genes_distribution.pdf"), width = 6, height = 6)
print(p)
dev.off()

# Heatmap ====
# Select all DEGs from deg_li, ordered by slope
# Select top 25 up-regulated and top 25 down-regulated genes based on slope
top_upregulated_genes <- deg_li %>%
  arrange(desc(slope)) %>%
  head(25) %>%
  pull(Gene)

top_downregulated_genes <- deg_li %>%
  arrange(slope) %>%
  head(25) %>%
  pull(Gene)

# Combine selected genes
selected_genes <- c(top_upregulated_genes, top_downregulated_genes)

# Subset count_matrix to include only these genes
heatmap_data <- count_matrix[rownames(count_matrix) %in% selected_genes, , drop = FALSE]

# Extract only SPF samples
spf_samples <- grep("SPF", colnames(heatmap_data), value = TRUE)

# Define the correct order of samples: SI → Ce → LI, each ordered by age (3, 17, 78)
ordered_samples <- c(
  grep("^SI_SPF_3", spf_samples, value = TRUE),
  grep("^SI_SPF_17", spf_samples, value = TRUE),
  grep("^SI_SPF_78", spf_samples, value = TRUE),
  grep("^Ce_SPF_3", spf_samples, value = TRUE),
  grep("^Ce_SPF_17", spf_samples, value = TRUE),
  grep("^Ce_SPF_78", spf_samples, value = TRUE),
  grep("^LI_SPF_3", spf_samples, value = TRUE),
  grep("^LI_SPF_17", spf_samples, value = TRUE),
  grep("^LI_SPF_78", spf_samples, value = TRUE)
)

# Ensure all samples exist in heatmap_data
ordered_samples <- ordered_samples[ordered_samples %in% colnames(heatmap_data)]

# Subset and reorder heatmap_data columns
heatmap_data <- heatmap_data[, ordered_samples]
heatmap_data <- heatmap_data[selected_genes, , drop = FALSE]

# Apply log2 transformation to normalize data
heatmap_data <- log2(heatmap_data + 1)  # Adding 1 to avoid log(0) issues

# Extract age and tissue from sample names for annotation
age_labels <- as.numeric(sub(".*_(\\d+)_\\d+$", "\\1", ordered_samples))
tissue_labels <- sub("_SPF_.*", "", ordered_samples)  # Extract tissue type

# Create annotation dataframe for columns
annotation_col <- data.frame(
  Age = factor(age_labels, levels = c(3, 17, 78)),  # Ensure correct ordering
  Tissue = factor(tissue_labels, levels = c("SI", "Ce", "LI"))  # Define tissue levels
)
rownames(annotation_col) <- ordered_samples

# Define annotation colors
ann_colors <- list(
  Age = c("3" = "lightblue", "17" = "orange", "78" = "red"),
  Tissue = c("SI" = "purple", "Ce" = "green", "LI" = "blue")
)

# Generate the heatmap with annotations
heatmap_plot <- pheatmap(
  heatmap_data,
  annotation_col = annotation_col,  # Add age and tissue annotations
  annotation_colors = ann_colors,   # Define colors for annotations
  cluster_rows = TRUE,              # Cluster genes
  cluster_cols = FALSE,             # Do not cluster samples (use custom order)
  scale = "row",                    # Normalize rows (genes)
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Top 25 Up- and Down-Regulated DEGs in LI Ordered by Slope",
  silent = TRUE
)

# Save the heatmap as a PDF
pdf(file = here::here(figure_dir, "Heatmap_top25_up_down_deg_li.pdf"), width = 8, height = 6)
grid.draw(heatmap_plot$gtable)  # Correct way to save pheatmap
dev.off()