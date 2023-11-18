# This repository hosts useful bash and R scripts for RNAseq and RRBSseq analysis.

## RRBSseq data processing
- [Trimgalore.sh](Trimgalore.sh): [reference](https://github.com/FelixKrueger/TrimGalore), score-trimming, adapter trimming. Different from typical trimming method, Trimgalore can automatically detect the existence of possible adapter sequences and remove them.
- [trimRRBSdiversityAdaptCustomers.py](trimRRBSdiversityAdaptCustomers.py): [reference](https://github.com/nugentechnologies/NuMetRRBS/blob/master/trimRRBSdiversityAdaptCustomers.py), trim adapter/indexes sequences from RRBSseq library.
- [trimRRBSdiversityAdaptCustomers.sh](./trimRRBSdiversityAdaptCustomers.sh): bash scripts for running trimRRBSdiversityAdaptCustomers.py in batch.
- [bismark.sh](./bismark.sh): [reference](https://github.com/FelixKrueger/Bismark),  map bisulfute-transformed sequence to genome.
- [bismark_methylation_extractor.sh](./bismark_methylation_extractor.sh): [reference](https://github.com/FelixKrueger/Bismark), extract methylation information.

## Identify differentially-methylated regions (DMR)
- [getMethylationCount.R](./getMethylationCount.R): Obtain the methylation counts (5hmc+5mc, 5hmc) from bismark extraction files.
- [get5mcCount.R](./get5mcCount.R): Obtain the real 5mc levels, ((5hmc + 5mc) - 5hmc).
- [getDMR.R](./getDMR.R): Obtain the differentially-methylated regons between germ-free (treatment) and specific-pathogen-free (control) condition.

## RNAseq data processing
- [salmon.sh](./salmon.sh): [reference](https://combine-lab.github.io/salmon/), map the RNAseq data to geme and determine the counts for each transcripts simultaneously. 