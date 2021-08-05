# Bacon

A comprehensive computational benchmarking framework for evaluating targeted chromatin conformation capture-specific methodologies

For more information, see our [Bacon](https://csuligroup.com/Bacon) webpage. 

## Contents
- [Prerequisites](#Prerequisites)
- [Example usage](#Example-usage)
- [Trouble shooting](#Trouble-shooting)

## Prerequisites

[Python](https://www.python.org/)(>=3.4.0)
[R](https://www.r-project.org/)(>=3.6.2)
[Homer](http://homer.ucsd.edu/homer/)
[bedtools](https://bedtools.readthedocs.io/en/latest/index.html)

## Example usage


### 1. UV Rate Calculation
Please prepare unfiltered alignment file of all PETs as input, the UV Rate will be printed on the screen.
```bash
# Run the commandline
UVRate_calculation.sh align_file.sam /tmpDir

```

### 2. PC Calculation
Please prepare at least one loop files (.bedpe format) and place them in a folder, high-quality ChIP-seq/CUT&Run peak file(.bed format).
```bash
# Run the commandline
Peak-occupancy.sh /directory_of_loop_folder peak.bed prefix /outDir

```

### 3. ES Calculation
Please prepare at least one loop files (.bedpe format) and place them in a folder, uniquely valid PET file (.bedpe format)
```bash
# Run the commandline
Enrichment_score.sh  /directory_of_loop_folder ValidPET.bedpe /outDir

```

### 4. ACC Calculation
Please prepare at least one loop files (.bedpe format) and place them in a folder, directory of gold standard(true) loop file, and directory of fale loops (you can find in folder /Gold_standard_loops).
Here we used dataset_1 as an example:
```bash
# Run the commandline
ACC_calculation.sh /directory_of_loop_folder /Gold_standard_loops/true /Gold_standard_loops/false dataset_1 /outDir

```

### 5. AR annotation
Please prepare at least one loop files (.bedpe format) and place them in a folder, ccres file of corresponding species (you can find in folder /AR_annotation), tss file of corresponding species (you can find in folder /AR_annotation), active histone mark ChIP-seq peaks (you can find in folder /AR_annotation), inactive histone mark ChIP-seq peaks (you can find in folder /AR_annotation).
```bash
# Run the commandline
AR_annotation.sh /directory_of_loop_folder /AR_annotation/hg19_ccres.bed /AR_annotation/GENCODEv19-TSSs.4k.bed /AR_annotation/peaks_k562/active /AR_annotation/peaks_k562/inactive prefix /outDir

```

## Trouble shooting
* Stdin error<br>
Please unzip all the annotation/loop files before running the script.

* Bedtools not found<br>
Please install bedtools and add the path into system.