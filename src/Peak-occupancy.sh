#!/bin/bash
#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# Calculate peak occupancy for peak-based methods
#===============
# bedtools need to be installed in system path (https://bedtools.readthedocs.io/en/latest/)
# loop = called loops in bedpe format
# peak = high quality ChIP-seq/CUT&RUN peak file
# prefix = the prefix of output file
# outDir = the directory of output files

loop_dir=$1
peak=$2
prefix=$3
outDir=$4

mkdir -p $outDir

if [ -f "$outDir/${prefix}_PC.txt" ]
then
    echo "Remove old file..."
    rm $outDir/${prefix}_PC.txt
fi

for loop in ${loop_dir}/*
do

file_name1="${loop##*/}" &&
file_name="${file_name1%.*}"

awk '{s1=$2-(10000-$3+$2)/2;e1=10000+s1;s2=$5-(10000-$6+$5)/2;e2=10000+s2;printf "%s\t%d\t%d\t%s\t%d\t%d\n",$1,s1,e1,$4,s2,e2}' $loop |  sed 's/-//' > $outDir/${file_name}.ext.loop.tmp
awk '{s1=$2-(10000-$3+$2)/2;e1=10000+s1;printf "%s\t%d\t%d\n",$1,s1,e1}' $peak |  sed 's/-//' > $outDir/${file_name}.ext.peak.tmp

awk '{print $1"\t"$2"\t"$3}' $outDir/${file_name}.ext.loop.tmp > $outDir/anchor1.tmp
awk '{print $4"\t"$5"\t"$6}' $outDir/${file_name}.ext.loop.tmp > $outDir/anchor2.tmp

cat $outDir/anchor1.tmp $outDir/anchor2.tmp | sort -u > $outDir/${file_name}.txt

bedtools intersect -c -a $outDir/${file_name}.ext.loop.tmp -b $outDir/${file_name}.ext.peak.tmp | awk '{if($4>0)print $1"\t"$2"\t"$3}' | sort -u > $outDir/${file_name}.overlap.txt
#bedtools intersect -c -a $peak -b $outDir/${file_name}.txt | awk '{print $1"\t"$2"\t"$3}' | sort -u > $outDir/${file_name}.overlap2.txt

anchor_num=$(awk 'END{print NR}' $outDir/${file_name}.txt)
overlap_anchor_count=$(< "$outDir/${file_name}.overlap.txt" wc -l)
overlap_read_count=$(awk '{sum+=$4} END{print sum}' $outDir/${file_name}.overlap.txt)

start_a=$(awk '{ total += $2 } END { print int(total/NR) }' $outDir/${file_name}.overlap.txt)
end_a=$(awk '{ total += $3 } END { print int(total/NR) }' $outDir/${file_name}.overlap.txt)
l_a=$(echo $end_a-$start_a | bc )

start_p=$(awk '{ total += $2 } END { print int(total/NR) }' $peak)
end_p=$(awk '{ total += $3 } END { print int(total/NR) }' $peak)
l_p=e$(echo $end_p-$start_p | bc)
 
ocRate=$(echo "scale=4; $overlap_anchor_count/$anchor_num" | bc)

pcRate=$(echo "scale=4; $ocRate/$l_a" | bc)

echo "peak occurancy for "${file_name}" : "$ocRate
echo "peak occurancy for "${file_name}" : "$ocRate >> $outDir/${prefix}_PC.txt

rm $outDir/${file_name}*
rm $outDir/*.tmp


done