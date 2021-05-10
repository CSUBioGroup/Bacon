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
ValidPet=$2
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


## Start calculating
## Convert loop/PETs files
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $InputLoop > $OutDir/InputLoop.tmp
echo "Convert input loop file successfully!"

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' $ValidPet > $OutDir/ValidPet.tmp
echo "Convert input valid PET file successfully!"

## Calculate density
Rscript Estimate_density.R $OutDir/InputLoop.tmp  $OutDir/ValidPet.tmp $OutDir
awk '{if($1==$4) sd=($3-$2)+($6-$5);s1=$2-sd;e1=$3+sd;s2=$5-sd;e2=$6+sd;print $1"\t"s1"\t"e1"\t"$4"\t"s2"\t"e2}' $InputLoop > $OutDir/InputLoop.tmp
loop_num=$(awk 'END{print NR}' $OutDir/InputLoop.tmp)
density_arry[0]="density"

N=20
for i in $(seq 1 $loop_num)
do
	chrom=$(awk -v var="$i" 'NR==var {print $1}' $OutDir/InputLoop.tmp)
	s1=$(awk -v var="$i" 'NR==var {print $2}' $OutDir/InputLoop.tmp)
	e1=$(awk -v var="$i" 'NR==var {print $3}' $OutDir/InputLoop.tmp)
	s2=$(awk -v var="$i" 'NR==var {print $5}' $OutDir/InputLoop.tmp)
	e2=$(awk -v var="$i" 'NR==var {print $6}' $OutDir/InputLoop.tmp)

	((i=i%N)); ((i++==0)) && wait
	awk -v chrom="$chrom" -v s1="$s1" -v e1="$e1" -v s2="$s2" -v e2="$e2" '{if($1==chrom && $4==chrom && $2>=s1 && $3<=e1 && $5>=s2 && $6<=e2)
	{density=density+1}};END{print chrom"\t"s1"\t"e1"\t"density}' $ValidPet >> $OutDir/InputLoop_density.tmp
	
	awk 'NR==FNR {a[$1]=$2;b[$1]=$3;c[$1$4]=$5;d[$1$4]=$6;next}
	{if($2>=a[$1$4] && $3<=b[$1$4] && $5>=c[$1$4] && $6<=d[$1$4])
	{density=density+1;print density}}'  $OutDir/InputLoop.tmp $ValidPet


done


printf "%s\n" "${density_arry[@]}" > $OutDir/Density_loop.txt
rm $OutDir/*.tmp



