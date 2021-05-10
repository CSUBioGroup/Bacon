#!/bin/bash
#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# function description
# calculate E-P percentage
# calculate AR
#===============

inputDir=$1
ccres=$2
tss=$3
active_peakDir=$4
inactive_peakDir=$5
outDir=$6

rm $outDir/result_*

for file in $inputDir/*
do

file_name1="${file##*/}" &&
file_name="${file_name1%.*}"

## Calculate E-P percentage

awk '{print $1"\t"$2"\t"$3}' $file > $outDir/${file_name}_loop1.tmp &&
awk '{print $4"\t"$5"\t"$6}' $file > $outDir/${file_name}_loop2.tmp &&

bedtools intersect -wa -c -a $outDir/${file_name}_loop1.tmp -b GENCODEv19-TSSs.4k.bed  > $outDir/${file_name}_loop1_promoter.tmp &&
bedtools intersect -wa -c -a $outDir/${file_name}_loop1.tmp -b K562_ccres.bed  > $outDir/${file_name}_loop1_enhancer.tmp &&

paste $outDir/${file_name}_loop1_promoter.tmp $outDir/${file_name}_loop1_enhancer.tmp > $outDir/${file_name}_loop1_merge.txt &&

awk '{if($4>$8){s="p"}else if($4<$8){s="e"}else if($4==$8 && $4==0){s="n"} else{s="e"};print $1"\t"$2"\t"$3"\t"s}' $outDir/${file_name}_loop1_merge.txt > $outDir/${file_name}_loop1_merge_a.txt &&

bedtools intersect -wa -c -a $outDir/${file_name}_loop2.tmp -b GENCODEv19-TSSs.4k.bed > $outDir/${file_name}_loop2_promoter.txt &&
bedtools intersect -wa -c -a $outDir/${file_name}_loop2.tmp -b K562_ccres.bed > $outDir/${file_name}_loop2_enhancer.txt &&

paste $outDir/${file_name}_loop2_promoter.txt $outDir/${file_name}_loop2_enhancer.txt > $outDir/${file_name}_loop2_merge.txt &&

awk '{if($4>$8){s="p"}else if($4<$8){s="e"}else if($4==$8 && $4==0){s="n"} else{s="e"};print $1"\t"$2"\t"$3"\t"s}' $outDir/${file_name}_loop2_merge.txt > $outDir/${file_name}_loop2_merge_a.txt &&

paste $outDir/${file_name}_loop1_merge_a.txt $outDir/${file_name}_loop2_merge_a.txt > $outDir/${file_name}_merge_all.txt &&

awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$4"_"$8}' $outDir/${file_name}_merge_all.txt > $outDir/${file_name}_merge_ep_anno.txt

awk '{if($7=="e_e" || $7=="e_p" || $7=="p_e" || $7=="p_p")print $0}' $outDir/${file_name}_merge_ep_anno.txt > $outDir/${file_name}_ep.txt

all_num=$(< "$outDir/${file_name}_merge_ep_anno.txt" wc -l) &&

rm $outDir/${file_name}_loop* &&
rm $outDir/${file_name}_merge* &&

e_e=$(awk '{if($7=="e_e")n=n+1}END{print n}' $outDir/${file_name}_ep.txt)
e_p=$(awk '{if($7=="e_p")n=n+1}END{print n}' $outDir/${file_name}_ep.txt)
p_p=$(awk '{if($7=="p_p")n=n+1}END{print n}' $outDir/${file_name}_ep.txt)
e_e_rate=$(echo "scale=4; $e_e/($all_num)" | bc)
e_p_rate=$(echo "scale=4; $e_p/($all_num)" | bc)
p_p_rate=$(echo "scale=4; $p_p/($all_num)" | bc)

echo -e $file_name"\t"$e_e_rate"\t"$e_p_rate"\t"$p_p_rate >> $outDir/result_percent.txt

echo -e $file_name": "$e_e_rate"; "$e_p_rate"; "$p_p_rate

## Calculate AA

awk '{print $1"\t"$2"\t"$3}' $outDir/${file_name}_ep.txt > $outDir/${file_name}_loop1.tmp &&
awk '{print $4"\t"$5"\t"$6}' $outDir/${file_name}_ep.txt > $outDir/${file_name}_loop2.tmp &&

for peak in $active_peakDir/*
do 
peak_name1="${peak##*/}" &&
peak_name="${peak_name1%.*}"
bedtools intersect -wa -c -a $outDir/${file_name}_loop1.tmp -b $peak  > $outDir/${file_name}_loop1_${peak_name}.tmp &&
bedtools intersect -wa -c -a $outDir/${file_name}_loop2.tmp -b $peak  > $outDir/${file_name}_loop2_${peak_name}.tmp 
done

paste $outDir/${file_name}_loop1_*.tmp > $outDir/${file_name}_loop1_merge.txt &&
paste $outDir/${file_name}_loop2_*.tmp > $outDir/${file_name}_loop2_merge.txt &&
awk '{if($4+$8+$12+$16>2){s="a"}else{s="o"};print $1"\t"$2"\t"$3"\t"s}' $outDir/${file_name}_loop1_merge.txt > $outDir/${file_name}_loop1_merge_a.txt &&
awk '{if($4+$8+$12+$16>2){s="a"}else{s="o"};print $1"\t"$2"\t"$3"\t"s}' $outDir/${file_name}_loop2_merge.txt > $outDir/${file_name}_loop2_merge_a.txt &&
paste $outDir/${file_name}_loop1_merge_a.txt $outDir/${file_name}_loop2_merge_a.txt > $outDir/${file_name}_merge_all.txt &&

awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$4"_"$8}' $outDir/${file_name}_merge_all.txt > $outDir/${file_name}_active.txt &&

rm $outDir/${file_name}_loop* &&
rm $outDir/${file_name}_merge* &&

awk '{print $1"\t"$2"\t"$3}' $outDir/${file_name}_ep.txt > $outDir/${file_name}_loop1.tmp &&
awk '{print $4"\t"$5"\t"$6}' $outDir/${file_name}_ep.txt > $outDir/${file_name}_loop2.tmp &&

for inact_peak in $inactive_peakDir/*
do
inpeak_name1="${inact_peak##*/}" &&
inpeak_name="${inpeak_name1%.*}"
bedtools intersect -wa -c -a $outDir/${file_name}_loop1.tmp -b $inact_peak  > $outDir/${file_name}_loop1_${inpeak_name}.txt &&
bedtools intersect -wa -c -a $outDir/${file_name}_loop2.tmp -b $inact_peak  > $outDir/${file_name}_loop2_${inpeak_name}.txt 
done

paste $outDir/${file_name}_loop1_*.txt > $outDir/${file_name}_loop1_merge.txt &&
paste $outDir/${file_name}_loop2_*.txt > $outDir/${file_name}_loop2_merge.txt &&
awk '{if($4>0){s="i"}else{s="o"};print $1"\t"$2"\t"$3"\t"s}' $outDir/${file_name}_loop1_merge.txt > $outDir/${file_name}_loop1_merge_a.txt &&
awk '{if($4>0){s="i"}else{s="o"};print $1"\t"$2"\t"$3"\t"s}' $outDir/${file_name}_loop2_merge.txt > $outDir/${file_name}_loop2_merge_a.txt &&
paste $outDir/${file_name}_loop1_merge_a.txt $outDir/${file_name}_loop2_merge_a.txt > $outDir/${file_name}_merge_all.txt &&
awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$6"\t"$7"\t"$4"_"$8}' $outDir/${file_name}_merge_all.txt > $outDir/${file_name}_inactive.txt &&

rm $outDir/${file_name}_loop* &&
rm $outDir/${file_name}_merge* &&

paste $outDir/${file_name}_active.txt $outDir/${file_name}_inactive.txt | awk '{if($7=="a_a"){s="a_a"}else if($7!="a_a" && $14=="i_i"){s="i_i"}else if(($7=="a_o" || $7=="o_a")&&($14=="o_o")){s="a_o"}else if(($7=="o_o")&&($14=="o_i" || $14=="i_o")){s="i_o"}else if($7==$14){s="o_o"}else if(($7=="o_a" || $7=="a_o")&&($14=="o_i" || $14=="i_o")){s="a_o"};print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"s}' \
 > $outDir/${file_name}_ac_inact.txt &&

rm $outDir/${file_name}_active* 
rm $outDir/${file_name}_inactive*

all_ac_num=$(< "$outDir/${file_name}_ac_inact.txt" wc -l) &&
all_file_num=$(< "$file" wc -l) &&

active_num=$(awk '{if($7=="a_a" ||$7=="a_o" ||$7=="o_a")n=n+1}END{print n}' $outDir/${file_name}_ac_inact.txt)
inactive_num=$(awk '{if($7=="i_i" ||$7=="i_o" ||$7=="o_i")n=n+1}END{print n}' $outDir/${file_name}_ac_inact.txt)
other_num=$(awk '{if($7=="o_o")n=n+1}END{print n}' $outDir/${file_name}_ac_inact.txt)
active_rate=$(echo "scale=4; $active_num/$all_ac_num" | bc)
inactive_rate=$(echo "scale=4; $inactive_num/$all_ac_num" | bc)
other_rate=$(echo "scale=4; $other_num/$all_ac_num" | bc)

echo -e $file_name"\t"$all_file_num"\t"$active_rate"\t"$inactive_rate"\t"$other_rate >> $outDir/result_AA.txt

echo $file_name": "$all_file_num";"$active_rate"; "$inactive_rate"; "$other_rate

done