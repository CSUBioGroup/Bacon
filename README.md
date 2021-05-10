# Bacon

A comprehensive computational benchmarking framework for evaluating targeted chromatin conformation capture-specific methodologies

For more information, see our [Bacon](https://csuligroup.com/Bacon) webpage. 

## Prerequisites

  [Python](https://www.python.org/)(>=3.4.0), 
  [R](https://www.r-project.org/)(>=3.6.2), 
  [juicer tool](https://github.com/aidenlab/juicer/wiki/Juicer-Tools-Quick-Start),
  [bedtools](https://bedtools.readthedocs.io/en/latest/index.html)

## Example usage


### 1. UV Rate Calculation
Please prepare unfiltered alignment file of all PETs as input, the UV Rate will be printed on the screen.
```bash
# Run the commandline
UVRate_calculation.sh align_file.sam /tmpDir

```

### 2. PC Calculation
Please prepare loop file (.bedpe format), high-quality ChIP-seq/CUT&Run peak file(.bed format).
```bash
# Run the commandline
Peak-occupancy.sh loop_file.bedpe peak.bed prefix /outDir

```

### 3. ES Calculation
Please prepare loop file (.bedpe format), uniquely valid PET file (.bedpe format) as inpu.
```bash
# Run the commandline
Enrichment_score.sh loop_file.bedpe ValidPET.bedpe /outDir

```

### 4. ACC Calculation
Please prepare loop file (.bedpe format), one significant gold standard loop file, and three negative gold loop files (check /Gold_standard_for_ACC).
Here we used K562-H3k27ac as an example:
```bash
# Run the commandline
cd /Gold_standard_for_ACC/K562-H3k27ac/count>5
ACC_calculation.sh loop_file.bedpe k562 /outDir

```
