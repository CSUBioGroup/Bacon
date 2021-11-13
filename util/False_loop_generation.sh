#!/bin/bash
#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# Function description
# False loops generation
#===============

genome_size=$1
gold=$2
hic_strong=$3
output=$4

if [ -f "$output" ]
then
    echo "Remove old file..."
    rm "$output"
fi

# number of loops for each chrom
count=$(< "$gold" wc -l)
num=$(echo "$(($count/10))")

# width of anchor, length of loop
start_g=$(awk '{ total += $2 } END { print int(total/NR) }' $gold)
end_g=$(awk '{ total += $3 } END { print int(total/NR) }' $gold)
start2_g=$(awk '{ total += $5 } END { print int(total/NR) }' $gold)

width=$(echo $end_g-$start_g | bc)
length=$(echo $start2_g-$end_g | bc)

echo "Generating total loops: $count"
echo "Generating loop width: $width"
echo "Generating loop length: $length"

while IFS=" " read f1 f2
do
for ((j=1;j<=10;j++))
do
echo "round: $j of $f1"
    for ((i=1;i<=$num/10;i++))
    do
    fold=$((10 ** $j))
    hgh=$f2
    RANGE1=$(echo $RANDOM % $hgh + $fold | bc)

    hgh2=$length
    RANGE2=$(echo $RANDOM % $hgh2 + 1 | bc)

    end1=$(echo $RANGE1 + $width | bc)
    start2=$(echo $RANGE1 + $width +$RANGE2 | bc)
    end2=$(echo $RANGE1+$width+$RANGE2+$width | bc)
    

    echo -e $f1"\t"$RANGE1"\t"$end1"\t"$f1"\t"$start2"\t"$end2 >> ${output}_tmp

    done
done
done <"$genome_size"

pairToPair -type neither -a ${output}_tmp -b $gold > ${output}_removegold_tmp.txt

pairToPair -type neither -a ${output}_removegold_tmp.txt -b $hic_strong > ${output}_output.txt

rm *.tmp

