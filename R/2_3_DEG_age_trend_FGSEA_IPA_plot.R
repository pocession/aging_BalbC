# Load necessary libraries ====
library(dplyr)
library(ggplot2)

# Dir ====
doc_dir <- here::here("Doc")
result_dir <- here::here("Results")
fgsea_dir <- here::here(result_dir, "Statistics", "FGSEA")
ipa_dir <- here::here(result_dir, "Statistics", "IPA")
figure_dir <- here::here(result_dir, "Figures")

# Read the FGSEA file
fgsea_results <- read.csv(here::here(fgsea_dir, "DEG_age_trend_li_fgsea_results.csv"))

# Step 4: Filter significant pathways (adjusted p-value < 0.05)
fgsea_results_filtered <- fgsea_results %>%
  filter(padj < 0.05) %>%
  arrange(padj)

# Step 1: Select the top 5 upregulated and top 5 downregulated pathways based on NES
top_upregulated <- fgsea_results_filtered %>%
  arrange(desc(NES)) %>%
  head(10)  # Top 5 most activated pathways

top_downregulated <- fgsea_results_filtered %>%
  arrange(NES) %>%
  head(10)  # Top 5 most suppressed pathways

# Combine the selected pathways
top_pathways <- bind_rows(top_upregulated, top_downregulated) %>%
  mutate(
    pathway_clean = gsub("^GOBP_", "", pathway),  # Remove "GO_BP_" prefix
    pathway_clean = gsub("_", " ", pathway_clean),  # Replace "_" with space
    neg_log10_padj = -log10(padj)  # Transform adjusted p-value
  )

# Step 2: Generate the dot plot with no inner grid lines
p <- ggplot(top_pathways, aes(x = NES, y = reorder(pathway_clean, NES), size = size, color = neg_log10_padj)) +
  geom_point(alpha = 0.8) +  # Use dots instead of bars
  scale_size(range = c(2, 6)) +  # Adjust dot size range (min, max)
  scale_color_gradient(low = "blue", high = "red") +  # Blue = less significant, Red = highly significant
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +  # Add padding to prevent touching edges
  labs(
    title = "Top Up- and Down-Regulated GO:BP Pathways (Top Genes in LI)",
    x = "Normalized Enrichment Score (NES)",
    y = "GO:BP Term",
    color = "-log10 Adjusted P-value",
    size = "Pathway Size"
  ) +
  theme_bw() +  # Apply theme_bw()
  theme(
    axis.text.y = element_text(size = 8),   # Reduce text size for readability
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# Save the plot as a PDF
pdf(file = here::here(figure_dir, "Dotplot_fgsea_top_up_down_deg_li_go_bp.pdf"), width = 8, height = 6)
print(p)
dev.off()

# IPA barplot ====
## TF ====
ipa_results <- read.delim(here::here(ipa_dir, "upstream_TF_deg_li.txt"), sep = "\t", skip = 1)

# Extract the number inside the parentheses from the "Mechanistic.Network" column
ipa_results <- ipa_results %>%
  dplyr::filter(Mechanistic.Network != " ") %>%
  mutate(Mechanistic.Network.Size = as.numeric(gsub(".*\\((\\d+)\\)", "\\1", Mechanistic.Network)))
  

# Step 1: Select the top 5 upregulated and top 5 downregulated master regulators based on Activation Z-score
top_upregulated <- ipa_results %>%
  arrange(desc(Activation.z.score)) %>%
  head(10)  # Top 5 most activated regulators

top_downregulated <- ipa_results %>%
  arrange(Activation.z.score) %>%
  head(10)  # Top 5 most suppressed regulators

# Combine the selected regulators
top_regulators <- bind_rows(top_upregulated, top_downregulated) %>%
  mutate(
    neg_log10_pval = -log10(p.value.of.overlap)  # Transform p-value of overlap
  )

# Step 2: Generate the dot plot using theme_bw()
p <- ggplot(top_regulators, aes(x = Activation.z.score, y = reorder(Upstream.Regulator, Activation.z.score), 
                                size = Mechanistic.Network.Size, color = neg_log10_pval)) +
  geom_point(alpha = 0.8) +  # Use dots instead of bars
  scale_size(range = c(2, 6)) +  # Adjust dot size range (min, max)
  scale_color_gradient(low = "blue", high = "red") +  # Blue = less significant, Red = highly significant
  labs(
    title = "Top 10 Up- and Down-Regulated Master Regulators from IPA",
    x = "Activation Z-score",
    y = "Master Regulator",
    color = "-log10 P-value of Overlap",
    size = "Mechanistic Network Size"
  ) +
  theme_bw() +  # Apply theme_bw()
  theme(
    axis.text.y = element_text(size = 8),   # Reduce text size for readability
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# Save the plot as a PDF
pdf(file = here::here(figure_dir, "Dotplot_ipa_top_up_down_TF.pdf"), width = 8, height = 6)
print(p)
dev.off()

## Cytokine ====
ipa_results <- read.delim(here::here(ipa_dir, "upstream_cytokine_deg_li.txt"), sep = "\t", skip = 1)

# Extract the number inside the parentheses from the "Mechanistic.Network" column
ipa_results <- ipa_results %>%
  dplyr::filter(Mechanistic.Network != " ") %>%
  mutate(Mechanistic.Network.Size = as.numeric(gsub(".*\\((\\d+)\\)", "\\1", Mechanistic.Network)))

# Step 1: Select the top 5 upregulated and top 5 downregulated master regulators based on Activation Z-score
top_upregulated <- ipa_results %>%
  arrange(desc(Activation.z.score)) %>%
  head(10)  # Top 5 most activated regulators

top_downregulated <- ipa_results %>%
  arrange(Activation.z.score) %>%
  head(10)  # Top 5 most suppressed regulators

# Combine the selected regulators
top_regulators <- bind_rows(top_upregulated, top_downregulated) %>%
  mutate(
    neg_log10_pval = -log10(p.value.of.overlap)  # Transform p-value of overlap
  )

# Step 2: Generate the dot plot using theme_bw()
p <- ggplot(top_regulators, aes(x = Activation.z.score, y = reorder(Upstream.Regulator, Activation.z.score), 
                                size = Mechanistic.Network.Size, color = neg_log10_pval)) +
  geom_point(alpha = 0.8) +  # Use dots instead of bars
  scale_size(range = c(2, 6)) +  # Adjust dot size range (min, max)
  scale_color_gradient(low = "blue", high = "red") +  # Blue = less significant, Red = highly significant
  labs(
    title = "Top 10 Up- and Down-Regulated Master Regulators from IPA",
    x = "Activation Z-score",
    y = "Master Regulator",
    color = "-log10 P-value of Overlap",
    size = "Mechanistic Network Size"
  ) +
  theme_bw() +  # Apply theme_bw()
  theme(
    axis.text.y = element_text(size = 8),   # Reduce text size for readability
    axis.text.x = element_text(size = 10),
    plot.title = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# Save the plot as a PDF
pdf(file = here::here(figure_dir, "Dotplot_ipa_top_up_down_cytokine.pdf"), width = 8, height = 6)
print(p)
dev.off()