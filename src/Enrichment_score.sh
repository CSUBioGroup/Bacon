#!/bin/bash

#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# Calculating Enrichment score for 
# cluster-based mthods
#======================
# Parameters for running
# InputLoop = loop file in bedpe format
# ValidPet = uniquely valid pet file in bedpe format
# Outdir = directory of output files
#======================

InputLoop=$1
Valid_dir=$2
Outdir=$3

#======================
# assign parameters
#======================

if [[ -d "$3" ]]; then
	rm -r $3;mkdir $3
	OutDir=$3
	echo "Output directory: $OutDir"
else
	OutDir=$3
	mkdir $OutDir
	echo "Output directory: $OutDir"
fi

for loop in ${InputLoop}/*
do

file_name1="${loop##*/}" &&
file_name="${file_name1%.*}"
IFS='_' read -ra ADDR <<< "$file_name"
name="${ADDR[0]}_${ADDR[1]}"
ValidPet="${Valid_dir}/${name}_valid.bedpe"

## Start calculating
## Convert loop/PETs files
awk '{if($1==$4)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $loop > $OutDir/InputLoop.tmp
echo "Convert input loop file successfully!"

awk '{if($1==$4)print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $ValidPet > $OutDir/ValidPet.tmp
echo "Convert input valid PET file successfully!"

## Calculate density
#Rscript Estimate_density.R $OutDir/InputLoop.tmp  $OutDir/ValidPet.tmp $OutDir
awk '{if($1==$4) sd=($3-$2)+($6-$5);s1=$2-sd;e1=$3+sd;s2=$5-sd;e2=$6+sd;print $1"\t"s1"\t"e1"\t"$4"\t"s2"\t"e2}' $loop > $OutDir/InputLoop.tmp
loop_num=$(awk 'END{print NR}' $OutDir/InputLoop.tmp)

pairToPair -type both -a $OutDir/InputLoop.tmp -b $OutDir/ValidPet.tmp > $OutDir/ValidPet_pair.tmp 

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $OutDir/ValidPet_pair.tmp | uniq -c >  $OutDir/${file_name}_loop_ES.txt

ES=$(awk '{ total += $1 } END { print total/NR }' $OutDir/${file_name}_loop_ES.txt)

echo $file_name" : "$ES >> $OutDir/global_ES_out.txt
echo $file_name" : "$ES 

#rm $OutDir/*.tmp

done

