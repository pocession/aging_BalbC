# This script is used to get methylation count from Bismark methylation extraction files

# Library ====
library(dplyr)
library(methylKit)
library(here)
library(stringr)

# Global variables ====
Year <- 2018

# Dir ====
sampleSheetDir <- here::here("Doc","sampleSheet")
inputDir <- here::here("Raw", "RRBS", paste0("cpg", "_", Year))
outputStatisticsDir <- here::here("Results", "Statistics")
outputTmpDir <- here::here("Results", "Tmp")
outputFiguresDir <- here::here("Results", "Figures")

# Sample sheet ====
sampleSheet <- read.csv(file.path(sampleSheetDir, paste0("RRBS",Year,".csv")))

## Append raw file names on the sample sheet
<<<<<<< HEAD
files <- data.frame(fname = list.files(inputDir), Sample = "")
for ( i in 1:nrow(files)) {
  files$Sample[i] <- str_split(files$fname[i],"_")[[1]][[1]]
}

sampleSheet <- sampleSheet |>
  dplyr::left_join(files, by = "Sample") |>
  dplyr::mutate(fpath = paste0(inputDir,"/",fname))
remove(files)
=======
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
>>>>>>> ef90439ef5fea4a42db2f8a7fd7a057145d9212d

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

# Get methylatio count ====
## Separate samples based on ages: 3, 17, 78
## Read methylation extraction files

sampleList <- sampleSheet |>
  dplyr::mutate(index = paste(Condition, Age_w, sep = "_"))

groupList <- sampleList |> dplyr::distinct(Age_w) |> unlist()

<<<<<<< HEAD
=======
## For 2018: treatment = rep(c(1,0), 4)
## For 2022: treatment = rep(c(1,0), 2)
>>>>>>> ef90439ef5fea4a42db2f8a7fd7a057145d9212d
myObjList <- list()
for (i in 1:length(groupList)) {
  sampleList_subset <- sampleList |> dplyr::filter(Age_w == groupList[[i]])
  myobj <- methylKit::methRead(as.list(sampleList_subset$fpath),sample.id = as.list(sampleList_subset$index), assembly="mm10",
<<<<<<< HEAD
                               treatment = rep(c(1,0), 2), mincov = 10, pipeline = "bismarkCytosineReport")
  myObjList[[groupList[i]]] <- myobj
}
remove(sampleList, i, myobj)


=======
                               treatment = rep(c(1,0), 4), mincov = 10, pipeline = "bismarkCytosineReport")
  myObjList[[groupList[i]]] <- myobj
}
remove(i, myobj)

## Filter the methylation reads
>>>>>>> ef90439ef5fea4a42db2f8a7fd7a057145d9212d
filteredMyObjList <- list()
for (i in 1:length(groupList)) {
  myobj <- myObjList[[groupList[i]]]
  filteredObj <- filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
  filteredMyObjList[[groupList[i]]] <- filteredObj
}
remove(i, myobj, filteredObj)

<<<<<<< HEAD
=======
## Summarize the methylation counts per 100bp
>>>>>>> ef90439ef5fea4a42db2f8a7fd7a057145d9212d
tileList <- list()
for (i in 1:length(groupList)) {
  filteredOobj <- filteredMyObjList[[groupList[i]]]
  tile <- tileMethylCounts(filteredOobj,win.size=100,step.size=100)
  meth=unite(tile, destrand=FALSE)
  tileList[[groupList[i]]] <- meth
}

# Save Meth data ====
for (i in 1:length(groupList)) {
  sampleList_subset <- sampleList |> dplyr::filter(Age_w == groupList[[i]])
  tissue <- sampleList_subset[i,]$Tissue
  write.csv(sampleList_subset,file.path(outputTmpDir, paste0("RRBS_",Year,"_", tissue,"_sampleSheet_", groupList[[i]],".csv")))
  tile.df <- tileList[[groupList[[i]]]]
  write.csv(tile.df, file.path(outputTmpDir, paste0("RRBS_", Year, "_", tissue, "_methTile100_", groupList[[i]],".csv")))
}

