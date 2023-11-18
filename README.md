# aging_BalbC
This project is to identify how murine genes in intestine are regulated by microbiota through methylation changes.

## Data sotrage
- All raw data and meta are stored in RIMLS server. The path is `/ceph/rimlsfnwi/data/cellbio/minoda/aging_BalbC`.
- All the downstream analysis, including DEA (differnetially-expression gene analysis), DMR (differentially-methylated regions), and microbiome 16S rRNA profileing, are temporarily saved in [Musa lab's google drive](https://drive.google.com/drive/folders/1WaWG5MWRS6cUywEomC_EXA32PuV3lVqf?usp=drive_link).

## RRBS
Reduced-representation bisulfite sequencing (RRBS-Seq ) experiments were conducted twice. The first experiment included small intestine samples and was conducted at 2018. The second experiment included large intestine sample and was conducted at 2022.

## RNANseqv1.1
poly-A mRNA sequencing experiments were conducted at 2018. The experiment included three tissues: SI, LI, CE, two conditions: SPF, GF, and three ages: 3w, 17w, 78w. The analysis has been conducted, please check this repository.
 
## microbiome
16SrDNA sequencing experiments were conducted at 2022-2023. The experiment included three tissues: SI, LI, CE, two conditions: SPF, GF, and three ages: 3w, 17w, 78w. Those files were not analyzed yet. 

## File structure
```
aging_BalbC
├── sampleSheet: sample sheet files
│   ├── *.xlsx
│ 
├── RRBS: all the methylation sequencing data
│   ├── fastq: raw fastq files
│   ├── ├── ├──*.fastq.gz
│   ├── ├── 2018
│   ├── ├── 2022
|   |
│   ├── trimmed: trimed fastq files
│   ├── ├── ├── 1. *.fq: trimed by Trimgalore
│   ├── ├── ├── 2. *.trimed.fq: trimed by a script provided by Nugen (This is the file used for bismark)
│   ├── ├── 2018
│   ├── ├── 2022
|   |
│   ├── bismapped : mapped files with bismark
│   ├── ├── ├── *_pe.bam: pair-end mapped files
│   ├── ├── ├──*_ambiguous_reads_*.fq.gz: sequences can't be mapped
│   ├── ├── ├── *.cytosine_context_summary.txt: this files are going to be used for downstream analysis. E.g Methylkit
│   ├── ├── 2018
│   ├── ├── 2022
|   |
|   ├── methylation_extraction: files containing methylated Cytosine under differnt contexts (CHH, Cpg, CHG)
│   ├── ├── 2018
│   ├── ├── 2022
|
├── RNAseqv1.1: all the RNA sequencing data
│   ├── fastq: raw fastq files
│   ├── ├── *.fastq.gz
|   |
│   ├── trimmed: trimed fastq files
│   ├── ├── *.fq.gz
|   |
│   ├── salmon: read counts 
│   ├── ├── sample
│   ├── ├── ├── salmon_quant
│   ├── ├── ├── ├── quant.sf : this file is going to ge used for downstream analysis. E.g DeSeq2
|
├── microbiome: all sequencing data of gut microbiomes
│   ├── fastq: raw fastq files
│   ├── ├── *.fastq.gz
```

## Change log
- 20231118: update DMR.
- 20231117: update 5hmc + 5mc, 5hmc, 5mc.