# Downstream analysis
1. Determination of 5-Methylcytosine (5mC) levels, 5-Hydroxymethylcytosine (5hmc) levels and transcripts per millions (TPM). 
2. The identification of differentially-expressed genes (DEG) and differentially-methylated regions (DMR). 
3. Integrated analysis (DEG+DMR).

# Analysis
Please read the following `.Rmd` file to know how we perform the downstream analysis.

## Unadjusted methylation counts
- [Unadjusted methylation counts of small intestine](./Doc/getUnadjustedMethylationCount2018.Rmd)
- [Unadjusted methylation counts of large intestine](./Doc/getUnadjustedMethylationCount2022.Rmd)

## Adjusted methylation counts
- [Adjusted methylation counts of small intestine](./Doc/getAdjustedMethylationCount2018.Rmd)
- [Adjusted methylation counts of large intestine](./Doc/getAdjustedMethylationCount2022.Rmd)

## Annotated differential methylated region (annotated DMR)

# Data storage
## Statistics data
In the [statistics repository](./Results/Statistics/), you can find the following analysis result.

- [5hmc + 5mc](./Results/Statistics/5hmc_5mc/): the 5mc + 5hmc levels for all samples.
- [5mc](./Results/Statistics/5mc/): 5mc levels for all samples.
- [5hmc](./Results/Statistics/5hmc/): 5hmc levels for all samples.
- [TMP](./Results/Statistics/TPM/): expression levles of genes for all samples, as shown in transcripts per millions (TPM). The files with prefix "vstTPM" are expression levels that have been applied by variance stabilizing transformation (vst). 
- [DEG](./Results/Statistics/DEG): differentially-expressed genes for all conditions.
- [DMR](./Results/Statistics/DMR/): differentially-methylated regions for all conditions, including 5hmc DMR and 5mc DMR.
- [annotated DMR](./Results/Statistics/DMR_annotated/): DMR annotated with the nearest genes as well as the distance to the nearest transcription starting site.

## Figures
In the [figures repository](./Results/Figures/), you can find the all figures.

## Tmp
In the [temporary repository](./Results/Tmp/), you can find all temporary files generated during the analysis. This is for developers only.
