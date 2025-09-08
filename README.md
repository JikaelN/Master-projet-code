# Qst–Fst Comparative Analysis under Neutral and Quantitative Trait Simulations

## Overview

This repository contains the full pipeline used to simulate, process, and analyze genetic data for comparing **Qst and Fst estimates** under different population structures and genetic architectures. The workflow combines:

- **Forward-time simulations (SLiM)**
- **Tree sequence processing (Python)**
- **Quantitative trait analysis (R)**

The repository provides:

- SLiM scripts for simulating *island* and *stepping-stone* models.  
- Python scripts to process simulated tree sequences into VCFs with added mutations and founder sampling.  
- R code implementing the **Guillaume–Whitlock (GW) framework** for Qst estimation and statistical testing.  
- Shell scripts for automated and parallelized execution on HPC clusters (SLURM).  

---

## Repository Structure


├── scripts/ # Core pipeline code

│ ├── Full_code_analysis_Qst_Fst.Rmd # Analysis and visualization of final results

│ ├── GW_Qst_Fst.r # Main algorithm (Guillaume–Whitlock framework)

│ ├── process_stepping_tree.py # Tree → VCF (stepping-stone model)

│ ├── process_tree_island.py # Tree → VCF (island model)

│ ├── run_batch_vcf.sh # Parallel VCF processing launcher

│ ├── run_one_vcf.sh # Single execution of tree → VCF

│ ├── run_one_seed.sh # Single SLiM simulation run

│ ├── run_slim_batch.sh # Array job: multiple SLiM simulations

│ ├── run_one_qst.sh # Single Qst–Fst analysis (R)

│ ├── run_qst_fst.sh # Batch Qst–Fst analysis (SLURM array)

│

├── slim_script/ # Forward simulations in SLiM

│ ├── island_neutral.slim # Island model simulation

│ ├── stepping_neutral.slim # Stepping-stone model simulation


I
---

## Installation

### Requirements

- **SLiM 4+** (forward-time population genetic simulator)  
- **Python 3.9+** with:  
  - `tskit`, `pyslim`, `msprime`, `numpy`, `pandas`  
- **R (≥ 4.2)** with:  
  - `tidyverse`, `hierfstat`, `JGTeach`, `vcfR`, `gaston`, `VGAM`, `boot`  


```bash
# Python environment
conda create -n qstfst python=3.10 numpy pandas pyslim tskit msprime
conda activate qstfst

# R environment (example)
module load R/4.3.1
Rscript -e 'install.packages(c("tidyverse","vcfR","boot","VGAM"))'
```

###Usage

1. Run Simulations


**Island model example:**
```bash
slim slim\_script/island\_neutral.slim
```

**Stepping-stone model:**
```bash
slim slim\_script/stepping\_neutral.slim
```

Or submit arrays on SLURM:
```bash
bash scripts/run\_slim\_batch.sh
```


2\. Process Tree Sequences → VCF

bash scripts/run\_batch\_vcf.sh



3\. Qst–Fst Analysis



Run a single replicate:



bash scripts/run\_one\_qst.sh





Or submit in parallel:



bash scripts/run\_qst\_fst.sh



4\. Visualization and Results



Open the R Markdown file for figure generation and summary statistics:



rmarkdown::render("scripts/Full\_code\_analysis\_Qst\_Fst.Rmd")



Authors \& Affiliations

* Author: Jikaël Ntoko, University of Lausanne (UNIL)
* Contributor (algorithmic component): Isabela Do’Ò (original algorithm adapted and extended)
