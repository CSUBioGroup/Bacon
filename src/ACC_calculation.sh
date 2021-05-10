#!/bin/bash
#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# function description
# Accuracy calculation
# bedtools need to be installed
#===============

input_loop=$1
prefix=$2
golden_sig=${prefix}"_golden_sig.txt"
golden_nosig1=${prefix}"_golden_nosig1.txt"
golden_nosig2=${prefix}"_golden_nosig2.txt"
golden_nosig3=${prefix}"_golden_nosig3.txt"
out_dir=$3

# extend golden sig loops to 5kb
rm $out_dir/results.txt

for i in $input_loop/*
do
file_name="${i##*/}" &&
input_num=$(< "$i" wc -l) &&
sig_num=$(< "$golden_sig" wc -l) &&

intra_num_1=$(awk '{if($1==$4)n=n+1}END{print n}' $i)
inter_num_1=$(awk '{if($1!=$4)n=n+1}END{print n}' $i)
intra_num_2=$(awk '{if($1==$4 && $7>2)n=n+1}END{print n}' $i)
inter_num_2=$(awk '{if($1!=$4 && $7>2)n=n+1}END{print n}' $i)
echo "intra_num_1: "$intra_num_1 >> $out_dir/results.txt
echo "inter_num_1: "$inter_num_1 >> $out_dir/results.txt
echo "intra_num_2: "$intra_num_2 >> $out_dir/results.txt
echo "inter_num_2: "$inter_num_2 >> $out_dir/results.txt

awk '{if($1==$4 && $7>3)print $0}' $i > $out_dir/intra_num2.tmp

awk '{s1=$2-(10000-$3+$2)/2;e1=10000+s1;s2=$5-(10000-$6+$5)/2;e2=10000+s2;printf "%s\t%d\t%d\t%s\t%d\t%d\n",$1,s1,e1,$4,s2,e2}' $golden_sig > $out_dir/$golden_sig.ext.tmp &&

# extend loop anchor to 10kb
awk '{s1=$2-(10000-$3+$2)/2;e1=10000+s1;s2=$5-(10000-$6+$5)/2;e2=10000+s2;if(s1<0)s1=0;if(s2<0)s2=0;printf "%s\t%d\t%d\t%s\t%d\t%d\n",$1,s1,e1,$4,s2,e2}' $out_dir/intra_num2.tmp > $out_dir/$prefix.ext.tmp &&
if [ "$input_num" -gt "$sig_num" ]
then
	sort -r -k 7 -n $out_dir/$prefix.ext.tmp > $out_dir/$prefix.ext.sort.tmp
	awk -v sig_num="${sig_num}" 'NR<=sig_num{print $0}' $out_dir/$prefix.ext.sort.tmp > $out_dir/$prefix1.ext.tmp 
else
	cp $out_dir/$prefix.ext.tmp $out_dir/$prefix1.ext.tmp 
fi

pairToPair -type either -a $out_dir/$golden_sig.ext.tmp -b $out_dir/$prefix1.ext.tmp | awk '{print $1"\t"$2"\t"$3}' | sort -u > $out_dir/$prefix.ext.sig.tmp &&
pairToPair -type either -a $golden_nosig1 -b $out_dir/$prefix1.ext.tmp | awk '{print $1"\t"$2"\t"$3}' | sort -u > $out_dir/$prefix.ext.nosig1.tmp &&
pairToPair -type either -a $golden_nosig2 -b $out_dir/$prefix1.ext.tmp | awk '{print $1"\t"$2"\t"$3}' | sort -u > $out_dir/$prefix.ext.nosig2.tmp &&
pairToPair -type either -a $golden_nosig3 -b $out_dir/$prefix1.ext.tmp | awk '{print $1"\t"$2"\t"$3}' | sort -u > $out_dir/$prefix.ext.nosig3.tmp &&

nosig_num=$(< "$golden_nosig1" wc -l) &&
TP=$(< "$out_dir/$prefix.ext.sig.tmp" wc -l) &&
FP=$(echo $sig_num-$TP | bc) &&
FN_round1=$(< "$out_dir/$prefix.ext.nosig1.tmp" wc -l)
TN_round1=$(echo $sig_num-$FN_round1 | bc)
FN_round2=$(< "$out_dir/$prefix.ext.nosig2.tmp" wc -l)
TN_round2=$(echo $sig_num-$FN_round2 | bc)
FN_round3=$(< "$out_dir/$prefix.ext.nosig3.tmp" wc -l)
TN_round3=$(echo $sig_num-$FN_round3 | bc)
FN_mean=$(((FN_round1+FN_round2+FN_round3)/3))
TN_mean=$(echo $sig_num-$FN_mean | bc)

PPV=$(echo "scale=4; $TP/($TP+$FP)" | bc)
TPR=$(echo "scale=4; $TP/($TP+$FN_mean)" | bc)
F_score=$(echo "scale=4; 2*$PPV*$TPR/($PPV+$TPR)" | bc)
tmp=$(echo "scale=4; ($FP+$FN_mean)*0.5" | bc)
F_score2=$(echo "scale=4; $TP/($TP+$tmp)" | bc)
Accuracy=$(echo "scale=4; ($TP+$TN_mean)/($TP+$FP+$TN_mean+$FN_mean)" | bc)

echo "TP: "$TP 
echo "sig_num: "$sig_num
echo "FP: "$FP
echo "FN_mean: "$FN_mean
echo "TN_mean: "$TN_mean

echo "PPV: "$PPV
echo "TPR: "$TPR
echo "tmp: "$tmp
#echo "F-score: "$F_score
echo "F-score: "$F_score2
echo "Accuracy: "$Accuracy

echo -e $file_name"\t"$TP"\t"$FP"\t"$FN_mean"\t"$TN_mean"\t"$PPV"\t"$TPR"\t"$F_score2"\t"$Accuracy >> $out_dir/results.txt &&

rm $out_dir/*.tmp 

done
wait

