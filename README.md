# aging_BalbC
This project is to identify how murine genes in intestine are regulated by microbiota through methylation changes.

## Data sotrage
All raw data and meta are stored in RIMLS server. The path is `/ceph/rimlsfnwi/data/cellbio/minoda/aging_BalbC`.

## RRBS
Reduced-representation bisulfite sequencing (RRBS-Seq ) experiments were conducted twice. The first experiment included small intestine samples and was conducted at 2018. The second experiment included large intestine sample and was conducted at 2022.

## File structure
```
aging_BalbC
├── meta files
├── RRBS: all the methylation sequencing data
│   ├── fastq: raw fastq files
│   ├── ├── 2018: *.fastq.gz
│   ├── ├── 2022: *.fastq.gz
│   ├── trimmed: trimed fastq files
│   ├── ├── *.fq: trimed by [Trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
│   ├── ├── *.trimed.fq: trimed by a script provided by [Nugen](https://github.com/nugentechnologies/NuMetRRBS/blob/master/README.md) (This is the file used for bismark).
│   ├── bismapped : mapped files with bismark
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
│   ├── ├── sample
│   ├── ├── ├── salmon_quant
│   ├── ├── ├── ├── quant.sf : this file is going to ge used for downstream analysis. E.g DeSeq2
├── microbome: all sequencing data of gut microbiomes
│   ├── fastq: raw fastq files
│   ├── ├── *.fastq.gz
```
