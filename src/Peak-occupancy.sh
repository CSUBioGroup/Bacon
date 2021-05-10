#!/bin/bash
# Calculate peak occupancy
# Created by LiTang
# tangli_csu@csu.edu.cn

# bedtools need to be installed in system path (https://bedtools.readthedocs.io/en/latest/)
# loop = called loops in bedpe format
# peak = high quality ChIP-seq/CUT&RUN peak file
# prefix = the prefix of output file
# outDir = the directory of output files

loop=$1
peak=$2
prefix=$3
outDir=$4

mkdir -p $outDir

awk '{print $1"\t"$2"\t"$3}' $loop > $outDir/anchor1.tmp
awk '{print $4"\t"$5"\t"$6}' $loop > $outDir/anchor2.tmp

cat $outDir/anchor1.tmp $outDir/anchor2.tmp | sort -u > $outDir/$prefix.txt

bedtools intersect -wo -a $outDir/$prefix.txt -b $peak | awk '{print $1"\t"$2"\t"$3}' | sort -u > $outDir/$prefix.overlap.txt

anchor_num=$(awk 'END{print NR}' $outDir/$prefix.txt)
overlap_anchor_num=$(awk 'END{print NR}' $outDir/$prefix.overlap.txt)
ocRate=$(echo "scale=4; $overlap_anchor_num/$anchor_num" | bc)

echo "peak occurancy for "$prefix" : "$ocRate

rm $outDir/*.tmp