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
gold_dir=$2
false_dir=$3
prefix=$4
golden_sig=${gold_dir}/${prefix}"_gold.bedpe"
false1=${false_dir}/${prefix}"_false1.bedpe"
false2=${false_dir}/${prefix}"_false2.bedpe"
false3=${false_dir}/${prefix}"_false3.bedpe"
out_dir=$5

if [ -f "$out_dir/${prefix}_results.txt" ]
then
    echo "Remove old file..."
    rm $out_dir/${prefix}_results.txt
fi

for i in $input_loop/*
do
file_name="${i##*/}" &&
input_num=$(< "$i" wc -l) &&
sig_num=$(< "$golden_sig" wc -l) &&

intra_num_1=$(awk -F "\t" '{if($1==$4)n=n+1}END{print n}' $i)
inter_num_1=$(awk '{if($1!=$4)n=n+1}END{print n}' $i)
intra_num_2=$(awk '{if($1==$4 && $7>2)n=n+1}END{print n}' $i)
inter_num_2=$(awk '{if($1!=$4 && $7>2)n=n+1}END{print n}' $i)

echo $file_name"; intra_num_: "$intra_num_1 >> $out_dir/${prefix}_results.txt

awk '{if($1==$4)print $0}' $i > $out_dir/intra_num2.tmp

awk '{s1=$2-(10000-$3+$2)/2;e1=10000+s1;s2=$5-(10000-$6+$5)/2;e2=10000+s2;printf "%s\t%d\t%d\t%s\t%d\t%d\n",$1,s1,e1,$4,s2,e2}' $golden_sig > $out_dir/$prefix.gold.ext.tmp &&

# extend loop anchor to 10kb
awk '{s1=$2-(10000-$3+$2)/2;e1=10000+s1;s2=$5-(10000-$6+$5)/2;e2=10000+s2;if(s1<0)s1=0;if(s2<0)s2=0;printf "%s\t%d\t%d\t%s\t%d\t%d\n",$1,s1,e1,$4,s2,e2}' $out_dir/intra_num2.tmp > $out_dir/$prefix.ext.tmp &&
#if [ "$input_num" -gt "$sig_num" ]
#then
#	sort -r -k 7 -n $out_dir/$prefix.ext.tmp > $out_dir/$prefix.ext.sort.tmp
#	awk -v sig_num="${sig_num}" 'NR<=sig_num{print $0}' $out_dir/$prefix.ext.sort.tmp > $out_dir/$prefix.1.ext.tmp 
#else
cp $out_dir/$prefix.ext.tmp $out_dir/$prefix.1.ext.tmp 
#fi

pairToPair -type both -a $out_dir/$prefix.gold.ext.tmp -b $out_dir/$prefix.1.ext.tmp | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' | sort -u > $out_dir/$prefix.ext.sig.tmp &&
pairToPair -type both -a $false1 -b $out_dir/$prefix.1.ext.tmp | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $out_dir/$prefix.ext.false1.tmp &&
pairToPair -type both -a $false2 -b $out_dir/$prefix.1.ext.tmp | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $out_dir/$prefix.ext.false2.tmp &&
pairToPair -type both -a $false3 -b $out_dir/$prefix.1.ext.tmp | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $out_dir/$prefix.ext.false3.tmp &&

false_num=sig_num
TP=$(< "$out_dir/$prefix.ext.sig.tmp" wc -l) &&
FN=$(echo $sig_num-$TP | bc) &&
FP_round1=$(< "$out_dir/$prefix.ext.false1.tmp" wc -l)
TN_round1=$(echo $sig_num-$FP_round1 | bc)
FP_round2=$(< "$out_dir/$prefix.ext.false2.tmp" wc -l)
TN_round2=$(echo $sig_num-$FP_round2 | bc)
FP_round3=$(< "$out_dir/$prefix.ext.false3.tmp" wc -l)
TN_round3=$(echo $sig_num-$FP_round3 | bc)
#FP_mean=$(((FP_round1+FP_round2+FP_round3)/3))
#TN_mean=$(echo $sig_num-$FP_mean | bc)
####################################################################################
PPV_1=$(echo "scale=4; $TP/($TP+$FP_round1)" | bc)
TPR_1=$(echo "scale=4; $TP/($TP+$FN)" | bc)
Accuracy_1=$(echo "scale=4; ($TP+$TN_round1)/($TP+$FP_round1+$TN_round1+$FN)" | bc)
####################################################################################
PPV_2=$(echo "scale=4; $TP/($TP+$FP_round2)" | bc)
TPR_2=$(echo "scale=4; $TP/($TP+$FN)" | bc)
Accuracy_2=$(echo "scale=4; ($TP+$TN_round2)/($TP+$FP_round2+$TN_round2+$FN)" | bc)
####################################################################################
PPV_3=$(echo "scale=4; $TP/($TP+$FP_round3)" | bc)
TPR_3=$(echo "scale=4; $TP/($TP+$FN)" | bc)
Accuracy_3=$(echo "scale=4; ($TP+$TN_round3)/($TP+$FP_round3+$TN_round3+$FN)" | bc)
####################################################################################

echo "FP_round1:"$FP_round1
echo "FP_round2:"$FP_round2
echo "FP_round3:"$FP_round3
echo "TN_round1:"$TN_round1
echo "TN_round2:"$TN_round2
echo "TN_round3:"$TN_round3
echo "PPV: "$PPV_1
echo "TPR: "$TPR_1
echo "Accuracy: "$Accuracy_1
echo "=============================="
echo "PPV: "$PPV_2
echo "TPR: "$TPR_2
echo "Accuracy: "$Accuracy_2
echo "=============================="
echo "PPV: "$PPV_3
echo "TPR: "$TPR_3
echo "Accuracy: "$Accuracy_3
echo "###################################################"

echo -e $TP"\t"$FP_round1"\t"$FN"\t"$TN_round1"\t"$PPV_1"\t"$TPR_1"\t"$Accuracy_1 >> $out_dir/${prefix}_results.txt &&
echo -e $TP"\t"$FP_round2"\t"$FN"\t"$TN_round2"\t"$PPV_2"\t"$TPR_2"\t"$Accuracy_2 >> $out_dir/${prefix}_results.txt &&
echo -e $TP"\t"$FP_round3"\t"$FN"\t"$TN_round3"\t"$PPV_3"\t"$TPR_3"\t"$Accuracy_3 >> $out_dir/${prefix}_results.txt &&

rm $out_dir/*.tmp 

done
wait

