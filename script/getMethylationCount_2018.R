# This script is used to get methylation count from Bismark methylation extraction files

# Library ====
library(dplyr)
library(methylKit)
library(here)
library(stringr)

# Global variables ====
YEAR <- 2018

# Dir ====
sampleSheetDir <- here::here("Doc","sampleSheet")
inputDir <- here::here("Raw", "RRBS", paste0("cpg", "_", YEAR))
outputStatisticsDir <- here::here("Results", "Statistics")
outputTmpDir <- here::here("Results", "Tmp")
outputFiguresDir <- here::here("Results", "Figures")

# Sample sheet ====
sampleSheet <- read.csv(file.path(sampleSheetDir, paste0("RRBS",YEAR,".csv")))

## Append raw file names on the sample sheet
## For 2018: files$Sample[i] <- str_split(files$fname[i],"-|_")[[1]][[2]]
## For 2022: files$Sample[i] <- str_split(files$fname[i],"_")[[1]][[1]]
files <- data.frame(fname = list.files(inputDir), Sample = "")
for ( i in 1:nrow(files)) {
  files$Sample[i] <- str_split(files$fname[i],"-|_")[[1]][[2]]
  files$Sample[i] <- files$Sample[i]
}

## For 2018
files <- files |>
  dplyr::mutate(ID = as.numeric(Sample)) |>
  dplyr::select(fname, ID)

## For 2018: dplyr::left_join(files, by = "ID") |>
## For 2022: dplyr::left_join(files, by = "Sample") |>
sampleSheet <- sampleSheet |>
  dplyr::left_join(files, by = "ID") |>
  dplyr::mutate(fpath = paste0(inputDir,"/",fname))
remove(files)

## Change the age
sampleSheet <- sampleSheet |>
  dplyr::mutate(Age_w = case_when(
    Age == "3 weeks" ~ "3",
    Age == "17 weeks" ~ "17",
    Age == "1.5 year" ~ "78"
  ))

## Get new sample index
sampleSheet <- sampleSheet |>
  dplyr::mutate(index = paste(Tissue, Condition, Age_w, OX, Replicate, sep = "_"))

# Get methylation count ====
## Read methylation extraction files
## GF = 1, SPF = 0
## For 2018: treatment = rep(c(1,1,0,0), 6)
## For 2022: treatment = rep(c(1, 0, 1, 0, 1, 0), 2)
## long process
methylObj <- methylKit::methRead(as.list(sampleSheet$fpath),sample.id = as.list(sampleSheet$index), 
                                 treatment = rep(c(1,1,0,0), 6), assembly="mm10", mincov = 10, pipeline = "bismarkCytosineReport")

## Filter the methylation count by coverage 
filteredObj <- filterByCoverage(methylObj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
remove(methylObj)

## Summarize the methylation counts per 100bp
## long process
## win.size = step.size; non-overlapping tiling
tile <- tileMethylCounts(filteredObj,win.size=100,step.size=100)

# meth <- unite(tile, destrand=FALSE)
remove(filteredObj)

# Save Meth data and sampleSheet ====
tissue <- sampleSheet[1,]$Tissue
saveRDS(tile, file.path(outputTmpDir, paste0("RRBS_", YEAR, "_", tissue, "_methTile100_", ".rds")))
write.csv(sampleSheet, file.path(outputTmpDir, paste0("RRBS_", YEAR, "_", tissue, "_sampleSheet_",".csv")))
# write.csv(tile, file.path(outputTmpDir, paste0("RRBS_", YEAR, "_", tissue, "_methTile100_",".csv")))

remove(tile)

