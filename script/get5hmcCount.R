# This script is used to get 5hmc count from methylation files that were already process by methylKit

library(dplyr)
library(methylKit)
library(here)
library(stringr)

# Global variables ====
# YEAR <- 2018
YEAR <- 2022

# TISSUE <- "SI"
TISSUE <- "LI"
group_list <- c("3", "17", "78")

# Dir ====
sampleSheetDir <- here::here("Doc","sampleSheet")
inputDir <- here::here("Results", "Tmp")
outputStatisticsDir <- here::here("Results", "Statistics")
outputTmpDir <- here::here("Results", "Tmp")
outputFiguresDir <- here::here("Results", "Figures")

# Sample sheet listt ====
sampleSheetList <- list()
for (i in 1:length(group_list)) {
  sample_sheet_subset <- read.csv(file.path(inputDir, paste0("RRBS","_", YEAR, "_", TISSUE, "_sampleSheet_", group_list[i], ".csv")))
  sample_sheet_subset <- sample_sheet_subset[,2:ncol(sample_sheet_subset)]
  sampleSheetList[[group_list[i]]] <- sample_sheet_subset
}
remove(i, sample_sheet_subset)

# Methylation call list ====
methList <- list()
for (i in 1:length(group_list)) {
  meth <- read.csv(file.path(inputDir, paste0("RRBS","_", YEAR, "_", TISSUE, "_methTile100_", group_list[i], ".csv")))
  meth <- meth[,2:ncol(meth)]
  methList[[group_list[i]]] <- meth
}
remove(i, meth)

# Get 5mc counts ====
## `OX+` = 5mc

for (i in 1:length(group_list)) {
  sample_sheet_subset <- sampleSheetList[[group_list[i]]] |>
    dplyr::filter(OX == "+")
  write.csv(sample_sheet_subset,file.path(outputTmpDir, paste0("RRBS_",YEAR,"_", TISSUE,"_5mc_sampleSheet_", group_list[[i]],".csv")))
}
remove(i, sample_sheet_subset)

# 2018: 5mc = 1, 3, 5, 7
# 2022: 5mc = 1, 2
for (i in 1:length(group_list)) {
  meth <- methList[[group_list[i]]] |>
    dplyr::select(chr, start, end, strand, 
                  coverage1, numCs1, numTs1,
                  coverage2, numCs2, numTs2
                  #coverage3, numCs3, numTs3,
                  #coverage5, numCs5, numTs5,
                  #coverage7, numCs7, numTs7
                  )
  write.csv(meth, file.path(outputTmpDir, paste0("RRBS_", YEAR, "_", TISSUE, "_5mc_methTile100_", group_list[[i]],".csv")))
}
remove(i, meth)
# Get 5hmc counts ====
## `OX-` = 5mc + 5hmc
## `OX-` - `OX+` = 5hmc
# 2018: 2-1, 4-3, 6-5, 8-7
# 2022: 3-1, 4-2
for (i in 1:length(group_list)) {
  meth <- methList[[group_list[i]]] |>
    dplyr::mutate(
      coverage3 = coverage3 - coverage1,
      numCs3 = numCs3 - numCs1,
      numTs3 = numTs3 - numTs1,
      coverage4 = coverage4 - coverage2,
      numCs4 = numCs4 - numCs2,
      numTs4 = numTs4 - numTs2
      #coverage2 = coverage2 - coverage1,
      #numCs2 = numCs2 - numCs1,
      #numTs2 = numTs2 - numTs1,
      #coverage4 = coverage4 - coverage3,
      #numCs4 = numCs4 - numCs3,
      #numTs4 = numTs4 - numTs3,
      #coverage6 = coverage6 - coverage5,
      #numCs6 = numCs6 - numCs5,
      #numTs6 = numTs6 - numTs5,
      #coverage8 = coverage8 - coverage7,
      #numCs8 = numCs8 - numCs7,
      #numTs8= numTs8 - numTs7
                  ) |>
    dplyr::select(chr, start, end, strand,
                  coverage3, numCs3, numTs3,
                  coverage4, numCs4, numTs4
                  #coverage4, numCs4, numTs4,
                  #coverage4, numCs4, numTs4,
                  #coverage2, numCs2, numTs2,
                  #coverage4, numCs4, numTs4,
                  #coverage6, numCs6, numTs6,
                  #coverage8, numCs8, numTs8
                  )
  write.csv(meth, file.path(outputTmpDir, paste0("RRBS_", YEAR, "_", TISSUE, "_5hmc_methTile100_", group_list[[i]],".csv")))
}
remove(i, meth)

