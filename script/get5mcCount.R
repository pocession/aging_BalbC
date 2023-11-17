# This script is used to get 5hmc and 5mc count from methylation files that were already process by methylKit

library(dplyr)
library(methylKit)
library(here)
library(stringr)

# Dir ====
sampleSheetDir <- here::here("Doc","sampleSheet")
inputTmpDir <- here::here("Results", "Tmp")
outputTmpDir <- here::here("Results", "Tmp")

# Global variables ====
# YEAR <- 2018
YEAR <- 2022

# TISSUE <- "SI"
TISSUE <- "LI"

# Sample sheet list ====
sampleSheet <- read.csv(file.path(inputTmpDir, paste0("RRBS_", YEAR, "_", TISSUE, "_sampleSheet_",".csv")))
sampleSheet  <- sampleSheet [,2:ncol(sampleSheet)]

sampleSheetTotal <- sampleSheet |>
  dplyr::filter(OX == "-") |>
  dplyr::mutate(fpath = file.path(outputTmpDir,paste0(paste0(index,"_total.csv"))))

sampleSheet5hmc <- sampleSheet |>
  dplyr::filter(OX == "+") |>
  dplyr::mutate(fpath = file.path(outputTmpDir,paste0(paste0(index,"_5hmc.csv"))))

sampleSheet5mc <- sampleSheet |>
  dplyr::filter(OX == "+") |>
  dplyr::mutate(fpath = file.path(outputTmpDir,paste0(paste0(index,"_5mc.csv"))))

# Read methylation call list (tile = 100) ====
tile <- readRDS(file.path(inputTmpDir, paste0("RRBS_", YEAR, "_", TISSUE, "_methTile100_", ".rds")))

# Confirm the sample ID and sampleSheet
sampleIdList <- getSampleID(tile)
for (i in 1:length(sampleIdList)) {
  if (sampleIdList[i] == sampleSheet[i,]$index)
    print("TRUE")
}

# Get 5hmc+5mc subset ====
# 2018: c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
# 2022: c(7, 8, 9, 10, 11, 12)
for (i in c(7, 8, 9, 10, 11, 12)) {
  sampleID <- sampleIdList[i]
  write.csv(tile[[i]], file.path(outputTmpDir,paste0(sampleID,"_total.csv")),row.names=FALSE,quote=FALSE)
}
remove(i, sampleID)

# Get 5hmc subset ====
# 2018: c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23)
# 2022: c(1, 2, 3, 4, 5, 6)
for (i in c(1, 2, 3, 4, 5, 6)) {
  sampleID <- sampleIdList[i]
  write.csv(tile[[i]], file.path(outputTmpDir,paste0(sampleID,"_5hmc.csv")),row.names=FALSE,quote=FALSE)
}
remove(i, sampleID)

# Get 5mc subset ====
# 2018: c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23)
# 2018: realMc <- methylKit::adjustMethylC(mc = tile[[i+1]], hmc = tile[[i]])
# 2022: c(1, 2, 3, 4, 5, 6)
# 2022: realMc <- methylKit::adjustMethylC(mc = tile[[i+6]], hmc = tile[[i]])
for (i in c(1, 2, 3, 4, 5, 6)) {
  sampleID <- sampleIdList[i]
  realMc <- methylKit::adjustMethylC(mc = tile[[i+6]], hmc = tile[[i]])
  write.csv(realMc, file.path(outputTmpDir,paste0(sampleID,"_5mc.csv")),row.names=FALSE,quote=FALSE)
}
remove(i, sampleID, realMc)



