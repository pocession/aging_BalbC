# aging_BalbC
This project is to identify how murine genes in intestine are regulated by microbiota through methylation changes.

## Data sotrage
All raw data and meta are stored in RIMLS server. The path is `/ceph/rimlsfnwi/data/cellbio/minoda/aging_BalbC`.

## File structure
```
aging_BalbC
├── meta files
├── RRBS1.4: all the methylation sequencing data
│   ├── fastq: raw fastq files
│   ├── *.fastq.gz
│   ├── RRBS_adapter_trimmed: trimed fastq files
│   ├── *.fq: trimed by [Trimgalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
│   ├── *.trimed.fq: trimed by a script provided by [Nugen](https://github.com/nugentechnologies/NuMetRRBS/blob/master/trimRRBSdiversityAdaptCustomers.py).
│   ├── bis_mapped : mapped files with [bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/)
│   ├── *_pe.bam: pair-end mapped files
│   ├── *_ambiguous_reads_*.fq.gz: sequences can't be mapped
│   ├── methylation_extraction: files containing methylated Cytosine under differnt contexts (CHH, Cpg, CHG)
│   ├── *.cytosine_context_summary.txt: this files are going to be used for downstream analysis. E.g.g [Methylkit](https://github.com/al2na/methylKit)
```
