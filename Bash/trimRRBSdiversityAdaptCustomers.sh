#!/bin/env bash
#SBATCH --partition=rimlsfnwi
#SBATCH --job-name=trimRRBSadapter
#SBATCH --output=./log/arr_%x-%A-%a.out
#SBATCH --error=./log/arr_%x-%A-%a.err
#SBATCH --time=24:00:00
#SBATCH --mem 100G
#SBATCH -c 16
#SBATCH --array=1-16

wd=/ceph/rimlsfnwi/data/cellbio/minoda/aging_BalbC
sub=RRBS

inputDir=$wd/$sub/trimmed
outputDir=$wd/$sub/trimmed

inputfile_list=($inputDir/*_R1_001.fastq.gz/*_R1_001_val_1.fq)
inputfile=${inputfile_list[$SLURM_ARRAY_TASK_ID-1]}

basename_temp=${inputfile%_R1_001_val_1.fq}
basename=${basename_temp##*/}

# S02140_S1_R1_001.fastq.gz

if [ -d "$outputDir/$basename" ]; then
        echo "outputDir/$basename exists."
        rm -r $outputDir/$basename
fi
mkdir $outputDir/$basename

# Read each raw fastq file 
INPUTFILENAME1=$inputfile
INPUTFILENAME2=${INPUTFILENAME1%_R1_001_val_1.fq}_R2_001_val_2.fq

# Trim diverse adapter sequence from 5' end
python $wd/$sub/trimRRBSdiversityAdaptCustomers.py -1 $INPUTFILENAME1 -2 $INPUTFILENAME2

OUTPUTFILENAME_temp1=${INPUTFILENAME1%.fq}_trimmed.fq
OUTPUTFILENAME_temp2=${INPUTFILENAME2%.fq}_trimmed.fq

fastqc $OUTPUTFILENAME_temp1
fastqc $OUTPUTFILENAME_temp2
