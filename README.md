# aging_BalbC
This project is to identify how murine genes in intestine are regulated by microbiota through methylation changes.

## Data sotrage
All raw data and meta are stored in RIMLS server. The path is `/ceph/rimlsfnwi/data/cellbio/minoda/aging_BalbC`.

## File structure
```
aging_BalbC
├── meta files
├── RRBSv1.4: all the methylation sequencing data
│   ├── fastq: raw fastq files
│   ├── ├── *.fastq.gz
│   ├── RRBS_adapter_trimmed: trimed fastq files
│   ├── ├── *.fq: trimed by [Trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
│   ├── ├── *.trimed.fq: trimed by a script provided by Nugen
│   ├── bis_mapped : mapped files with bismark
│   ├── ├── *_pe.bam: pair-end mapped files
│   ├── ├── *_ambiguous_reads_*.fq.gz: sequences can't be mapped
│   ├── methylation_extraction: files containing methylated Cytosine under differnt contexts (CHH, Cpg, CHG)
│   ├── ├── *.cytosine_context_summary.txt: this files are going to be used for downstream analysis. E.g Methylkit
├── RNAseqv1.1: all the methylation sequencing data
│   ├── fastq: raw fastq files
│   ├── ├── *.fastq.gz
│   ├── trimmed: trimed fastq files
│   ├── ├── *.fq.gz
│   ├── salmon: read counts 
│   ├── sample
│   ├── ├── salmon_quant
│   ├── ├── ├── quant.sf : this file is going to ge used for downstream analysis. E.g DeSeq2
├── microbome: all sequencing data of gut microbiomes
│   ├── fastq: raw fastq files
│   ├── ├── *.fastq.gz
```
