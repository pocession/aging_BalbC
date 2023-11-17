# get DMR
# This script is used to calculate differentially-methylated regions
## read meth first
methylObj <- methylKit::methRead(as.list(sampleSheetTotal$fpath),sample.id = as.list(sampleSheetTotal$index), 
                                 treatment = rep(c(1,0), 6), sep = ",", assembly="mm10")