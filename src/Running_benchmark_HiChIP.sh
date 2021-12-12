#!/bin/bash
#===============
# author: 
# Li Tang
# Central South University
# tangli_csu@csu.edu.cn
#===============
# Function description
# Benchmark of HiChIP tools
#===============

# Running FitHiChIP tool
input_dir=$1

for data in $input_dir/*
do

bash /home/tangli/tools/FitHiChIP/FitHiChIP_HiCPro.sh -C /home/tangli/tools/FitHiChIP/${data} &

done

# Running HiChipper

hichipper --out ${input_dir}/hichipper_out_C_S  ${input_dir}/example_COMBINED_SELF.yaml
hichipper --out ${input_dir}/hichipper_out_E_S  ${input_dir}/example_EACH_SELF.yaml
hichipper --out ${input_dir}/hichipper_out_C_A  ${input_dir}/example_COMBINED_ALL.yaml
hichipper --out ${input_dir}/hichipper_out_C_S  ${input_dir}/example_COMBINED_SELF.yaml
hichipper --out ${input_dir}/hichipper_out_peak  ${input_dir}/example_peak.yaml



