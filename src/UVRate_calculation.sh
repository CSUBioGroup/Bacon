#!/bin/bash

#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# Calculating UV Rate for datasets
#======================
# Parameters for running
# align_file = input the alignment file of PETs (.sam)
# tmpDir = directory to store tmp files


align_file=$1
tmpDir=$2

# Picard remove duplication
picardCMD="java -jar /tools/picard.jar"


## Extract file name
file_name1="${align##*/}" &&
prefix="${file_name1%.*}"

## Sort by coordinate
$picardCMD SortSam I=$align_file O=$tmpDir/$prefix.sorted.sam SORT_ORDER=coordinate
## mark duplicates
$picardCMD MarkDuplicates I=$tmpDir/$prefix.sorted.sam O=$tmpDir/$prefix.sorted.dupMarked.sam METRICS_FILE=$tmpDir/$prefix.picard.dupMark.txt
## remove duplicates
$picardCMD MarkDuplicates I=$tmpDir/$prefix.sorted.sam O=$tmpDir/$prefix.sorted.rmDup.sam REMOVE_DUPLICATES=true METRICS_FILE=$tmpDir/$prefix.picard.rmDup.txt
## Filter and keep the mapped read pairs
samtools view -S -F 0x04 $tmpDir/$prefix.sorted.rmDup.sam > $tmpDir/$prefix.sorted.rmDup.mapped.sam

rmdup_bedpe=$tmpDir/$prefix.sorted.rmDup.mapped.sam

all_count=$(awk 'END{print NR}' $align_file)
rmdup_count=$(awk 'END{print NR}' $rmdup_bedpe)

UVRate=$(echo "scale=4; $rmdup_count/$all_count" | bc)

echo "UVRate:"${UVRate}