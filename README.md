# aging_BalbC
This project is to identify how murine genes in intestine are regulated by microbiota through methylation changes.

## Data sotrage and analysis
- All raw data and meta files are stored in RIMLS server. The path is `/ceph/rimlsfnwi/data/cellbio/minoda/aging_BalbC`. Please see [readmeRawData.md](./readmeRawData.md) for more information.
- The pre-processed raw data are all stored in the [pre-processed repository](./Results/Preprocessed/). Please see [readmePreprocessedData.md](./readmePreprocessData.md) for more information.
- All the downstream analysis, including DEA (differnetially-expression gene analysis, from bulk RNAseq data), DMR (differentially-methylated regions, from bulk RRBSseq data), and microbiome 16S rRNA profileing, are saved in [the result repository](./Results/). Please see the [readmeDownstreamData.md](./readmeDownstreamData.md) for more information.

## Experiments
### RNANseq
poly-A mRNA sequencing experiments were conducted at 2018. The experiment included three tissues: SI, LI, CE, two conditions: SPF, GF, and three ages: 3w, 17w, 78w.

### RRBS
Two reduced-representation bisulfite sequencing (RRBS-Seq) experiments were conducted. The first experiment included small intestine samples and was conducted at 2018. The second experiment included large intestine sample and was conducted at 2022.
 
### microbiome
16SrDNA sequencing experiments were conducted at 2022-2023. The experiment included three tissues: SI, LI, CE, two conditions: SPF, GF, and three ages: 3w, 17w, 78w. Those files were not analyzed yet. 

## Change log
- 20240509: update functions. Also update annotated DMR files to [DMR_annotated](./Results/Statistics/DMR_annotated/) repository. Codes of this repository are also going to public. All data and presentation files are not shown.
- 20240225: re-organize the readme files and scripts.
- 20240224: upload annotated DMR files to [DMR_annotated](./Results/Statistics/DMR_annotated/) repository.
- 20240218: upload all raw files and temporary results. Re-organized the [Results](./Results/) repository.
- 20231118: update DMR.
- 20231117: update 5hmc + 5mc, 5hmc, 5mc.