# Raw data
All the raw data are saved in `/ceph/rimlsfnwi/data/cellbio/minoda/aging_BalbC`. See below for more information.

# File structure
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