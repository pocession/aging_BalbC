#!/bin/env bash
#SBATCH --partition=rimlsfnwi
#SBATCH --job-name=Trimgalore
#SBATCH --output=./log/arr_%x-%A-%a.out
#SBATCH --error=./log/arr_%x-%A-%a.err
#SBATCH --time=24:00:00
#SBATCH --mem 100G
#SBATCH -c 16
#SBATCH --array=1-16

wd=/ceph/rimlsfnwi/data/cellbio/minoda/aging_BalbC
sub=RRBS

inputDir=$wd/$sub/fastq/2022
outputDir=$wd/$sub/trimmed

inputfile_list=($inputDir/*_R1_001.fastq.gz)
inputfile=${inputfile_list[$SLURM_ARRAY_TASK_ID-1]}

basename_temp=${inputfile%_R1.fastq.gz}
basename=${basename_temp##*/}

# S02140_S1_R1_001.fastq.gz
if [ -d "$outputDir/$basename" ]; then
        echo "outputDir/$basename exists."
        rm -r $outputDir/$basename
fi
mkdir $outputDir/$basename

# Read each raw fastq file 
INPUTFILENAME1=$inputfile
INPUTFILENAME2=${INPUTFILENAME1%_R1_001.fastq.gz}_R2_001.fastq.gz

# Bismark can only recognize fastq, the output file should be .fastq
trim_galore --fastqc --rrbs --dont_gzip --paired --cores 4 \
-o $outputDir/$basename \
$INPUTFILENAME1 $INPUTFILENAME2
