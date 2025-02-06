# Perform FGSEA for the age trend genes in LI

# Load necessary libraries ====
library(dplyr)
library(ggplot2)

# Dir ====
doc_dir <- here::here("Doc")
gmt_dir <- here::here(doc_dir, "gmt")
result_dir <- here::here("Results")
fgsea_dir <- here::here(result_dir, "Statistics", "FGSEA")
deg_age_trend_dir <- here::here(result_dir, "Statistics", "DEG_age_trend")
ipa_dir <- here::here(result_dir, "Statistics", "IPA")
figure_dir <- here::here(result_dir, "Figures")

# Read the nature aging gene ====
deg_li <- read.csv(here::here(deg_age_trend_dir, "deg_li_age_trend.csv"))

# FGSEA bar plot ====
# Step 1: Select the top 50 DEGs from deg_li, ordered by slope
ranked_gene_list <- deg_li %>%
  arrange(desc(slope)) %>%
  dplyr::select(Gene, slope) %>%
  deframe()  # Convert to named numeric vector (gene names as names)


# Step 2: Choose pathways
# Load the GO:BP pathways from a .gmt file
pathways <- fgsea::gmtPathways(here::here(gmt_dir, "m5.go.bp.v2024.1.Mm.symbols.gmt"))

# Step 4: Perform fgsea analysis
fgsea_results <- fgsea(
  pathways = pathways, 
  stats = ranked_gene_list, 
  minSize = 15,  # Minimum genes in a pathway
  maxSize = 500, # Maximum genes in a pathway
  nperm = 10000  # Number of permutations for significance
)

# Convert list columns to character format for saving
fgsea_results_clean <- fgsea_results %>%
  mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = "; ")))

# Save the FGSEA file
write.csv(fgsea_results_clean, here::here(fgsea_dir, "DEG_age_trend_li_fgsea_results.csv"))