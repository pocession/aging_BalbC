#!/bin/env bash
#SBATCH --partition=rimlsfnwi
#SBATCH --job-name=bismark
#SBATCH --output=./log/arr_%x-%A-%a.out
#SBATCH --error=./log/arr_%x-%A-%a.err
#SBATCH --time=60:00:00
#SBATCH --mem 100G
#SBATCH -c 16
#SBATCH --array=1-24

# 2018, 24; 2022, 16

wd=/ceph/rimlsfnwi/data/cellbio/minoda/aging_BalbC
sub=RRBS
# year="2022"
year="2018"

inputDir=$wd/$sub/bismapped
outputDir=$wd/$sub/methylation_extraction/$year

genome_dir="/ceph/rimlsfnwi/data/cellbio/mhlanga/thsieh/mm10"

# 2018: RRBS-14_1_val_1_trimmed_bismark_bt2_pe.bam
# 2022: S02144_S5_R1_001_val_1_trimmed_bismark_bt2_pe.bam

# Read the raw fastq file list
# 2018
INPUTFILES_temp=($inputDir/$year/RRBS-*_1_val_1_trimmed_bismark_bt2_pe.bam)

# 2022
# INPUTFILES_temp=($inputDir/S021*/*_R1_001_val_1_trimmed_bismark_bt2_pe.bam)

# Read each raw fastq file 
INPUTFILES="${INPUTFILES_temp[$SLURM_ARRAY_TASK_ID-1]}"
INPUTFILENAME1=${INPUTFILES}

# Generate basename
# 2018
basename_temp=${INPUTFILENAME1%_1_val_1_trimmed_bismark_bt2_pe.bam}

# 2022
# basename_temp=${INPUTFILENAME1%_R1_001_val_1_trimmed_bismark_bt2_pe.bam}

basename=${basename_temp##*/}

if [ -d "$outputDir/$basename" ]; then
        echo "outputDir/$basename exists."
        rm -r $outputDir/$basename
fi
mkdir $outputDir/$basename

# Extract methyaiton
bismark_methylation_extractor --gzip --bedGraph $INPUTFILENAME1 -o $outputDir/$basename
bismark_methylation_extractor --gzip --comprehensive $INPUTFILENAME1 -o $outputDir/$basename
bismark_methylation_extractor --gzip --bedGraph --cytosine_report --genome_folder $genome_dir $INPUTFILENAME1 -o $outputDir/$basename
