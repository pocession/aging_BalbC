# aging_BalbC
This project is to identify how murine genes in intestine are regulated by microbiota through methylation changes.

## Data sotrage
- All raw data and meta files are stored in RIMLS server. The path is `/ceph/rimlsfnwi/data/cellbio/minoda/aging_BalbC`. Please see [README_rawdata.md](./README_rawdata.md) for more information.
- The processed raw data are all stored in the [processed repository](./Processed/). Please see [README_processeddata.md](./Processed/README_processeddata.md) for more information.
- All the downstream analysis, including DEA (differnetially-expression gene analysis), DMR (differentially-methylated regions), and microbiome 16S rRNA profileing, are saved in [the result repository](./Results/). Please see the [README_downstream.md](./Results/README_downstream.md) for more information.

## Experiments and data analysis
### RRBS
Reduced-representation bisulfite sequencing (RRBS-Seq ) experiments were conducted twice. The first experiment included small intestine samples and was conducted at 2018. The second experiment included large intestine sample and was conducted at 2022.

### RNANseqv1.1
poly-A mRNA sequencing experiments were conducted at 2018. The experiment included three tissues: SI, LI, CE, two conditions: SPF, GF, and three ages: 3w, 17w, 78w.
 
### microbiome
16SrDNA sequencing experiments were conducted at 2022-2023. The experiment included three tissues: SI, LI, CE, two conditions: SPF, GF, and three ages: 3w, 17w, 78w. Those files were not analyzed yet. 

## Change log
- 20240218: upload all raw files and temporary results. Re-organized the [Results](./Results/) repository.
- 20231118: update DMR.
- 20231117: update 5hmc + 5mc, 5hmc, 5mc.