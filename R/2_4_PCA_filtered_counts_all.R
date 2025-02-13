# Load necessary libraries ====
library(dplyr)
library(ggplot2)
library(ggfortify)
library(ggrepel)
source(here::here("R", "Functions.R"))

# Dir ====
doc_dir <- here::here("Doc")
result_dir <- here::here("Results")
deg_dir <- here::here(result_dir, "Statistics", "DEG")
count_dir <- here::here(result_dir, "Statistics", "RNA_count")
figure_dir <- here::here(result_dir, "Figures")

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

# Count ====
count <- read.csv(here::here(result_dir, "Statistics", "RNA_count", "vst_count_all_organ_all_age.csv"))

# Remove duplicate rows based on the SYMBOL column
count_unique <- count %>%
  distinct(SYMBOL, .keep_all = TRUE)

rownames(count_unique) <- count_unique$SYMBOL
count_unique <- count_unique[,2:ncol(count_unique)]

# Plot ====
df <- count_unique
pca <- prcomp(t(df))
PC <- c(1,2)
p <- autoplot(pca, data = sampleSheet, x = PC[1], y = PC[2], 
              colour = "Weeks.factor",
              shape = "Organ_GF.SPF", size = 4, label = T, label.size = 0,
              label.repel=F, 
              main = "PCA all") + 
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_shape_manual(values=c(0,15,1,16,2,17)) + scale_colour_manual(values = c("red", "blue", "green"))
p
ggsave(here::here(figure_dir, "PCA_all.pdf"),plot=p)
dev.off()

