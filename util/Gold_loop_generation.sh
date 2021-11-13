#!/bin/bash
#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# Function description
# gold standard loops generation
#===============

eqtl=$1
crispr=$2
gencode_v19=$3
candidate=$4
out_dir=$5

# extend genomic loci of variants
awk '{s1=$2-(5000-$3+$2)/2;e1=5000+s1;printf "%s\t%d\t%d\n",$1,s1,e1}' $eqtl > $out_dir/eqtl_ext.tmp &&

# extract genomic loci of Gencode gene
awk 'NR==FNR{a[$4]=$1"\t"$2"\t"$3} $4 in a{print $1"\t"$2"\t"$3"\t"a[$4]}' $out_dir/eqtl_ext.tmp > $out_dir/eqtl_loop.txt &&

cat $out_dir/eqtl_loop.txt $crispr > $out_dir/eqtl_crispr_loop.txt

awk '{print $1"\t"$2"\t"$3}' $out_dir/eqtl_crispr_loop.txt > $out_dir/eqtl_crispr_loop_a1.tmp
awk '{print $4"\t"$5"\t"$6}' $out_dir/eqtl_crispr_loop.txt > $out_dir/eqtl_crispr_loop_a2.tmp

awk '{print $1"\t"$2"\t"$3}' $candidate > $out_dir/$candidate_a1.tmp
awk '{print $4"\t"$5"\t"$6}' $candidate > $out_dir/$candidate_a2.tmp

cat $out_dir/eqtl_crispr_loop_a1.tmp $out_dir/eqtl_crispr_loop_a2.tmp > $out_dir/eqtl_crispr_loop_a.tmp

bedtools intersect -wa -c -a $out_dir/$candidate_a1.tmp -b $out_dir/eqtl_crispr_loop_a1.tmp | awk '{print $NF}' > $out_dir/$candidate_a1_a1.tmp
bedtools intersect -wa -c -a $out_dir/$candidate_a1.tmp -b $out_dir/eqtl_crispr_loop_a2.tmp | awk '{print $NF}' > $out_dir/$candidate_a1_a2.tmp
bedtools intersect -wa -c -a $out_dir/$candidate_a2.tmp -b $out_dir/eqtl_crispr_loop_a1.tmp | awk '{print $NF}' > $out_dir/$candidate_a2_a1.tmp
bedtools intersect -wa -c -a $out_dir/$candidate_a2.tmp -b $out_dir/eqtl_crispr_loop_a2.tmp | awk '{print $NF}' > $out_dir/$candidate_a2_a2.tmp

paste $out_dir/$candidate_a1_a1.tmp $out_dir/$candidate_a1_a2.tmp $out_dir/$candidate_a2_a1.tmp $out_dir/$candidate_a2_a2.tmp > $out_dir/$candidate_list.tmp

# define function to calculate p-value 
function statcom ( mq, mi, mj, mb )
{
    zz = 1
    mz = zz
    mk = mi
    while ( mk <= mj ) {
        zz = zz * mq * mk / ( mk - mb)
        mz = mz + zz
        mk = mk + 2
    }
    return mz
}

function studpval ( mt , mn )
{
    PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
    if ( mt < 0 )
        mt = -mt
    mw = mt / sqrt(mn)
    th = atan2(mw, 1)
    if ( mn == 1 )
        return 1.0 - th / (PI/2.0)
    sth = sin(th)
    cth = cos(th)
    if ( mn % 2 == 1 )
        return 1.0 - (th+sth*cth*statcom(cth*cth, 2, mn-3, -1))/(PI/2.0)
    else
        return 1.0 - sth * statcom(cth*cth, 1, mn-3, -1)
}


awk '{print system("studpval" ($1,$2,$3,$4),4)}' $out_dir/$candidate_list.tmp > $out_dir/$candidate_pval.tmp
paste $candidate $out_dir/$candidate_pval.tmp > $out_dir/$candidate_pval.txt

awk '{if($NF<0.05 && $(NF-1)<0.05)print $0}' $out_dir/$candidate_pval.txt > $out_dir/$candidat_gold.txt