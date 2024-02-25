#!/bin/env bash
#SBATCH --partition=rimlsfnwi
#SBATCH --job-name=bismark
#SBATCH --output=./log/arr_%x-%A-%a.out
#SBATCH --error=./log/arr_%x-%A-%a.err
#SBATCH --time=60:00:00
#SBATCH --mem 120G
#SBATCH -c 16
#SBATCH --array=1-16

wd=/ceph/rimlsfnwi/data/cellbio/minoda/aging_BalbC
sub=RRBS

inputDir=$wd/$sub/trimmed
outputDir=$wd/$sub/bismapped

inputfile_list=($inputDir/*_R1_001.fastq.gz/*_R1_001_val_1_trimmed.fq)
inputfile=${inputfile_list[$SLURM_ARRAY_TASK_ID-1]}

basename_temp=${inputfile%_R1_001_val_1_trimmed.fq}
basename=${basename_temp##*/}

# S02140_S1_R1_001.fastq.gz

if [ -d "$outputDir/$basename" ]; then
        echo "outputDir/$basename exists."
        rm -r $outputDir/$basename
fi
mkdir $outputDir/$basename

# Read the raw fastq file list
genome_dir="/ceph/rimlsfnwi/data/cellbio/mhlanga/thsieh/mm10/bismark"

# Read each raw fastq file 
INPUTFILENAME1=${inputfile}
INPUTFILENAME2=${INPUTFILENAME1%_R1_001_val_1_trimmed.fq}_R2_001_val_2_trimmed.fq

# Map to bisulfite converted genome
# Perform paired end mapping
bismark --genome $genome_dir -o $outputDir/$basename --un --ambiguous -1 $INPUTFILENAME1 -2 $INPUTFILENAME2
