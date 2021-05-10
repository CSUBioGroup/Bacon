#!/bin/bash
# Generate APA for significant loops
# Created by LiTang
# tangli_csu@csu.edu.cn

# Process bedpe file
# Parameters for running
# bedpe = input bedpe file with at least 7 columns: chr1 str1 end1 chr2 str2 end2 pval/count
# hic = the path of .hic file for the respective genome
# lift = [1,0] 1 run liftover; 0 don't run liftover
# prefix = the prefix of output files and temp files
# script = the script directory should include scripts: liftOverBedpe.py, liftOver
# juicer = the path of juicer tools with the format of .jar 
# chain = the chain file for liftover
# outdir = the directory of output files

bedpe=$1
hic=$2
lift=$3
prefix=$4
juicer=$5
outdir=$6

#################### preprocess bedpe file #############################################
echo -e "preprocess bedpe file..."

awk '{printf "%s\t%d\t%d\t%s\t%d\t%d\t%.9f\n", $1, $2, $3, $4, $5, $6, $7}' $bedpe  > $outdir/$prefix.temp.bedpe

judge_count=$(awk 'NR==1{if($7>1)print 1;else print 0}' $outdir/$prefix.temp.bedpe)

if [[ $judge_count -eq 1 ]]; then
	echo "preprocess bedpe file based on loop counts..."
	sort -nrk 7 $outdir/$prefix.temp.bedpe |  awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | uniq -u > $outdir/$prefix.temp.sort.bedpe
else
	echo "preprocess bedpe file based on p-value..."
	sort -nk 7 $outdir/$prefix.temp.bedpe |  awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' | uniq -u > $outdir/$prefix.temp.sort.bedpe
fi

awk '{n=n+1;print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\tpair_"n}' $outdir/$prefix.temp.sort.bedpe > $outdir/$prefix.temp.sort.id.bedpe
head -n 10000 $outdir/$prefix.temp.sort.id.bedpe > $outdir/$prefix.processed.bedpe

if [[ $lift -eq 1 ]]; then
	script=$7
	chain=$8

#################### liftover mm10 to mm9 ###############################################
echo -e "liftover mm10 to mm9..."

python $script/liftOverBedpe.py --lift $script/liftOver --chain $chain --i $outdir/$prefix.temp.sort.id.bedpe --o $outdir/$prefix.temp.mm10.bedpe
cp $outdir/$prefix.temp.mm10.bedpe $outdir/$prefix.processed.bedpe

fi

rm $outdir/*.temp.*
################ Running juicer tools to generate APA ###################################
java -Xmx2g -jar $juicer apa $hic $outdir/$prefix.processed.bedpe $outdir/$prefix

