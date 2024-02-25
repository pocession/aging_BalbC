# get DMR
# This script is used to calculate differentially-methylated regions

# library ====
library(dplyr)
library(methylKit)
library(here)
library(stringr)

# Dir ====
sampleSheetDir <- here::here("Results", "Tmp")
inputTmpDir <- here::here("Results", "Tmp")
outputTmpDir <- here::here("Results", "Tmp")

# Global variables ====
# YEAR <- 2018
YEAR <- 2022

# TISSUE <- "SI"
TISSUE <- "LI"

group_list <- c(3, 17, 78)

# Get sample sheet ====
sampleSheet <- read.csv(file.path(sampleSheetDir, paste0("RRBS_", YEAR, "_", TISSUE, "_sampleSheet_",".csv")))
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


## read meth ====
## 5hmc ====
# Create sample sheet list
sampleSheetDf <- sampleSheet5hmc
sampleSheetList <- list()
for (i in 1:length(group_list)) {
  sampleSheet_subset <- sampleSheetDf |>
    dplyr::filter(Age_w == group_list[i])
  index <- as.character(group_list[i])
  sampleSheetList[[index]] <- sampleSheet_subset
}
remove(i, index,sampleSheetDf, sampleSheet_subset)

## Read data
## 2018: treatment = c(1, 0, 1, 0)
## 2022: treatment = c(1, 0)
methylObjList <- list()
for (i in 1:length(group_list)) {
  methylObj <- methylKit::methRead(as.list(sampleSheetList[[i]]$fpath),sample.id = as.list(sampleSheetList[[i]]$index), 
                                   treatment =  c(1, 0), sep = ",", 
                                   pipeline=list(fraction=FALSE,chr.col=1,start.col=2,end.col=3,
                                                 strand.col=4, coverage.col=5,freqC.col=8), 
                                   assembly="mm10")
  index <- as.character(group_list[i])
  methylObjList[[index]] <- methylObj
}
remove(i, index, methylObj)

## Get DMR ====
dmrObjList <- list()
for (i in 1:length(group_list)) {
  index <- as.character(group_list[[i]])
  meth <- unite(methylObjList[[index]], destrand = FALSE)
  dmrObjList[[index]] <- calculateDiffMeth(meth)
}
remove(i, index, meth)

## Save into files ====
for (i in 1:length(group_list)) {
  index <- as.character(group_list[[i]])
  write.csv(dmrObjList[[index]], file.path(outputTmpDir, paste0("DMR_5hmc_", YEAR, "_", TISSUE, "_", index, ".csv")), row.names = FALSE, quote = FALSE)
}
remove(i, index)

## 5mc ====
# Create sample sheet list
sampleSheetDf <- sampleSheet5mc
sampleSheetList <- list()
for (i in 1:length(group_list)) {
  sampleSheet_subset <- sampleSheetDf |>
    dplyr::filter(Age_w == group_list[i])
  index <- as.character(group_list[i])
  sampleSheetList[[index]] <- sampleSheet_subset
}
remove(i, index,sampleSheetDf, sampleSheet_subset)

## Read data
## 2018: treatment = c(1, 0, 1, 0)
## 2022: treatment = c(1, 0)
methylObjList <- list()
for (i in 1:length(group_list)) {
  methylObj <- methylKit::methRead(as.list(sampleSheetList[[i]]$fpath),sample.id = as.list(sampleSheetList[[i]]$index), 
                                   treatment =  c(1, 0), sep = ",", 
                                   pipeline=list(fraction=FALSE,chr.col=1,start.col=2,end.col=3,
                                                 strand.col=4, coverage.col=5,freqC.col=8), 
                                   assembly="mm10")
  index <- as.character(group_list[i])
  methylObjList[[index]] <- methylObj
}
remove(i, index, methylObj)

## Get DMR ====
dmrObjList <- list()
for (i in 1:length(group_list)) {
  index <- as.character(group_list[[i]])
  meth <- unite(methylObjList[[index]], destrand = FALSE)
  dmrObjList[[index]] <- calculateDiffMeth(meth)
}
remove(i, index, meth)

## Save into files ====
for (i in 1:length(group_list)) {
  index <- as.character(group_list[[i]])
  write.csv(dmrObjList[[index]], file.path(outputTmpDir, paste0("DMR_5mc_", YEAR, "_", TISSUE, "_", index, ".csv")), row.names = FALSE, quote = FALSE)
}
remove(i, index)


