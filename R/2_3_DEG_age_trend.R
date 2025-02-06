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


# Function to detect genes with significant age-dependent changes and extract slopes
detect_age_trend <- function(count_data) {
  # Extract age from column names (assumed format: Tissue_Treatment_Age_Repeat)
  ages <- as.numeric(sub(".*_(\\d+)_\\d+$", "\\1", colnames(count_data)))
  
  # Initialize results dataframe
  results <- data.frame(Gene = rownames(count_data), p_value = NA, slope = NA, stringsAsFactors = FALSE)
  
  # Loop through each gene and perform regression analysis with log10-transformed counts
  for (i in 1:nrow(count_data)) {
    gene_expression <- as.numeric(count_data[i, ])
    
    # Apply log10 transformation (adding 1 to avoid log(0) issues)
    log_expression <- log2(gene_expression + 1)
    
    # Fit linear regression model
    model <- lm(log_expression ~ ages)
    summary_model <- summary(model)
    
    # Store p-value and slope of age variable
    results$p_value[i] <- coef(summary_model)[2, 4]
    results$slope[i] <- coef(summary_model)[2, 1]  # Slope of age
  }
  
  # Adjust p-values for multiple testing
  results$padj <- p.adjust(results$p_value, method = "BH")
  
  # Return significant genes
  results <- results %>%
    filter(padj < 0.05) %>%
    arrange(padj)
  
  return(results)
}

# Apply the function to each tissue category
deg_si <- detect_age_trend(count_si)
deg_li <- detect_age_trend(count_li)
deg_ce <- detect_age_trend(count_ce)

# Save the results
write.csv(deg_si, here::here(deg_age_trend_dir, "deg_si_age_trend.csv"), row.names = FALSE)
write.csv(deg_li, here::here(deg_age_trend_dir, "deg_li_age_trend.csv"), row.names = FALSE)
write.csv(deg_ce, here::here(deg_age_trend_dir, "deg_ce_age_trend.csv"), row.names = FALSE)