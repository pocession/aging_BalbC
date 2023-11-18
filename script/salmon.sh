#!/bin/bash
#$ -t 1-72
#$ -o /osc-fs_home/hsieh/RNAseqv1.1/log
#$ -e /osc-fs_home/hsieh/RNAseqv1.1/log
#

# set path
export PATH="/osc-fs_home/hsieh/salmon-1.6.0_linux_x86_64/bin:$PATH"
# set wd dir
wd=/osc-fs_home/hsieh/RNAseqv1.1

# Read the dir list
INPUTFILES_temp=($wd/trimmed/*trimmed.fq.gz)

# salmon index
index="/osc-fs_home/hsieh/genome/salmon_partial_sa_index"

# Assign the inputdir and inputfilename 
inputfile="${INPUTFILES_temp[$SGE_TASK_ID - 1]}"
basename_temp=${inputfile%_trimmed.fq.gz}
basename=${basename_temp##*/}

cd $wd
cd salmon

if [ -d "$basename" ]; then
        echo "$basename exists"
	rm  -r $basename
fi

mkdir $basename
cd $basename

if [ -d "salmon_quant" ]; then
	rm  -r salmon_quant
fi

# salmon count
salmon quant -i $index -l A -r $inputfile --validateMappings -o salmon_quant
