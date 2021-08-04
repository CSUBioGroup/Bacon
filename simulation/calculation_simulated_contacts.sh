#!/bin/bash
#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# function description
# Calculate ChIP coverage for simulated contacts
#===============

ChIP_peak=$1
Simu_contact=$2
cc_cutoff=$3
prefix=$4
out_dir=$5


awk '{print $1"\t"$2"\t"$3}' $Simu_contact > $out_dir/Simu.a1.tmp
awk '{print $4"\t"$5"\t"$6}' $Simu_contact > $out_dir/Simu.a2.tmp
#awk -v '{print $0}' $Simu_contact > $out_dir/Simu.$prefix.txt

cat $out_dir/Simu.a1.tmp $out_dir/Simu.a2.tmp | sort -u > $out_dir/Simu.a.tmp

bedtools intersect -wa -a $out_dir/Simu.a.tmp -b $ChIP_peak | sort -u > $out_dir/Simu.cc.tmp

count_simu=$(< "$out_dir/Simu.a.tmp" wc -l)
count_cc=$(< "$out_dir/Simu.cc.tmp" wc -l)
Rate_cc=$(echo "scale=3; $count_cc*100/$count_simu" | bc)

echo "count_simu: "$count_simu
echo "count_cc: "$count_cc
echo "Rate_cc: "$Rate_cc
