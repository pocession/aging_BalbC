# The pre-processed data
The raw data of RNAseq and RRBSseq were further pre-processed with the following bash and Python scripts. The pre-processed data are stored in the [Preprocessed data repository](./Results/Preprocessed/).

# Preprocess
## RRBSseq data processing
- [Trimgalore.sh](./Bash/Trimgalore.sh): [reference](https://github.com/FelixKrueger/TrimGalore), score-trimming, adapter trimming. Different from typical trimming method, Trimgalore can automatically detect the existence of possible adapter sequences and remove them.
- [trimRRBSdiversityAdaptCustomers.py](./Python/trimRRBSdiversityAdaptCustomers.py): [reference](https://github.com/nugentechnologies/NuMetRRBS/blob/master/trimRRBSdiversityAdaptCustomers.py), trim adapter/indexes sequences from RRBSseq library.
- [trimRRBSdiversityAdaptCustomers.sh](./Bash/trimRRBSdiversityAdaptCustomers.sh): bash scripts for running trimRRBSdiversityAdaptCustomers.py in batch.
- [bismark.sh](./Bash/bismark.sh): [reference](https://github.com/FelixKrueger/Bismark),  map bisulfute-transformed sequence to genome.
- [bismark_methylation_extractor.sh](./Bash/bismark_methylation_extractor.sh): [reference](https://github.com/FelixKrueger/Bismark), extract methylation information.

## RNAseq data processing
- [salmon.sh](./Bash/salmon.sh): [reference](https://combine-lab.github.io/salmon/), map the RNAseq data to geme and determine the counts for each transcripts simultaneously. 

# Data strorage
## RRBSseq
- [RRBS](./Results/Preprocessed/RRBS/): This repository hosts the CpG methylation calls, which are processed by bismark. The data of small intestine are stored in the [small intestine repository](./Results/Preprocessed/RRBS/cpg_2018/). The data of large intestine is tore in the [large intestine repository](./Results/Preprocessed/RRBS/cpg_2022/).

## RNAseq
- [RNAseq](./Results/Preprocessed/RNAseq/): This repository hosts the raw counts of transcripts, which are processed by salmon.