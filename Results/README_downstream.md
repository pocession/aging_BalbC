# Downstream analysis
The repository hosts analysis, including 5-Methylcytosine (5mC) levels, 5-Hydroxymethylcytosine (5hmc) levels, transcripts per millions (TPM), differentially-expressed genes (DEG), differentially-methylated regions (DMR), and integrated analysis (DEG+DMR).

## Statistics
Host all data frames.

- [Total](./Statistics/Total/): the 5mc + 5hmc levels for all samples.
- [5mc](./Statistics/5mc/): 5mc levels for all samples.
- [5hmc](./Statistics/5hmc/): 5hmc levels for all samples.
- [TMP](./Statistics/TPM/): expression levles of genes for all samples, as shown in transcripts per millions (TPM). The files with prefix "vstTPM" are expression levels that have been applied by variance stabilizing transformation (vst). 
- [DEG](./Statistics/DEG): differentially-expressed genes for all conditions.
- [DMR](./Statistics/DMR/): differentially-methylated regions for all conditions, including 5hmc DMR and 5mc DMR.

## Figures
Host all figures.

## Tmp
Host temporary files generated during the analysis. For developers only.
