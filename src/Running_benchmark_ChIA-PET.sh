#!/bin/bash
#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# Function description
# Benchmark of ChIA-PET tools
#===============

raw_dir=$1
prefix=$2
out_dir=$3


# Set your own configurations

config=$4
source ${config}
echo ${thread}
 
# Running CPT2
mkdir ${out_dir}/CPT2
ChIA-PET2 -g ${bwa_idx} -b ${chrom_hg19} -f ${raw_dir}/*_1.fastq.gz -r ${raw_dir}/*_2.fastq.gz -A ${linkerA} -B ${linkerB} -o ${out_dir}"/CPT2" -n ${prefix} -t ${thread} -Q 10 -m 0 -d ${short} -e 3 -l 18 &
P1=$!

# Running ChIAPoP
mkdir ${out_dir}/ChIAPoP
Rscript --vanilla ChIAPoP_running.R ${raw_dir} ${linkerA} ${linkerB} ${out_dir}/ChIAPoP &
P2=$!

# Running mango
mkdir ${out_dir}/mango
Rscript /home/tangli/tools/mango/mango/mango.R --fastq1 ${raw_dir}/*_1.fastq --fastq2 ${raw_dir}/*_2.fastq --prefix ${prefix} --outdir ${out_dir}"/mango"  --linkerA ${linkerA} --linkerB ${linkerB} --argsfile /home/tangli/tools/mango/mango/argfile.txt --reportallpairs TRUE --maxlength 100 --minlength 10 --gsize hs --chromexclude chrM --stages 3:5 &
P3=$!

# Running CPTv.3
mkdir ${out_dir}/CPT3
if [[ -f "linker.tmp" ]]; then
    rm linker.tmp
fi
echo ${linkerA} >> linker.tmp
echo ${linkerB} >> linker.tmp
java -jar /home/tangli/tools/ChIA-PET_Tool_V3/ChIA-PET.jar --mode 0 --fastq1 ${raw_dir}/*_1.fastq --fastq2 ${raw_dir}/*_2.fastq \
--linker linker.tmp --minimum_linker_alignment_score 5 --GENOME_INDEX ${bwa_idx} --GENOME_LENGTH ${GENOME_LENGTH} \
--CHROM_SIZE_INFO ${CHROM_SIZE_INFO} --CYTOBAND_DATA ${CYTOBAND_DATA} --SPECIES 1 --output ${out_dir}/CPT3 \
--prefix ${prefix} --MAPPING_CUTOFF 5 --thread 6 --minimum_tag_length 10 --minSecondBestScoreDiff 1 &
--MIN_COVERAGE_FOR_PEAK 2
P4=$!

# Running CID
mkdir ${out_dir}/CID
java -Xmx32G -jar gem.jar CID --data ${out_dir}/CPT2/${prefix}.rmdup.bedpe --g ${chrom_hg19} --micc 2 --out ${out_dir}/CID/${prefix}

# Running cLoops
mkdir ${out_dir}/cLoops
cLoops -f ${out_dir}/CPT2/${prefix}.rmdup.bedpe -o ${out_dir}/cLoops/${prefix} -s -m 1

wait
