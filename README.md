Qst–Fst Comparative Analysis under Neutral and Quantitative Trait Simulations

Overview



This repository contains the full pipeline used to simulate, process, and analyze genetic data for comparing Qst and Fst estimates under different population structures and genetic architectures. The workflow combines forward-time simulations (SLiM), tree sequence processing (Python), and quantitative trait analysis (R).



The repository provides:



SLiM scripts for simulating island and stepping-stone models.



Python scripts to process simulated tree sequences into VCFs with added mutations and founder sampling.



R code implementing the Guillaume–Whitlock (GW) framework for Qst estimation and statistical testing.



Shell scripts for automated and parallelized execution on HPC clusters (SLURM).



Repository Structure

├── scripts/                          # Core pipeline code

│   ├── Full\_code\_analysis\_Qst\_Fst.Rmd   # Analysis and visualization of final results

│   ├── GW\_Qst\_Fst.r                     # Main algorithm (Guillaume–Whitlock framework)

│   ├── process\_stepping\_tree.py         # Tree → VCF (stepping-stone model)

│   ├── process\_tree\_island.py           # Tree → VCF (island model)

│   ├── run\_batch\_vcf.sh                 # Parallel VCF processing launcher

│   ├── run\_one\_vcf.sh                   # Single execution of tree → VCF

│   ├── run\_one\_seed.sh                  # Single SLiM simulation run

│   ├── run\_slim\_batch.sh                # Array job: multiple SLiM simulations

│   ├── run\_one\_qst.sh                   # Single Qst–Fst analysis (R)

│   ├── run\_qst\_fst.sh                   # Batch Qst–Fst analysis (SLURM array)

│

├── slim\_script/                      # Forward simulations in SLiM

│   ├── island\_neutral.slim              # Island model simulation

│   ├── stepping\_neutral.slim            # Stepping-stone model simulation



Installation

Requirements



* SLiM 4+ (forward-time population genetic simulator)



* Python 3.9+ with: tskit, pyslim, msprime, numpy, pandas



* R (≥ 4.2) with: tidyverse, hierfstat, JGTeach, vcfR, gaston, VGAM, boot



Usage

1. Run Simulations

\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_\_

Island model example:



slim slim\_script/island\_neutral.slim





Stepping-stone model:



slim slim\_script/stepping\_neutral.slim





Or submit arrays on SLURM:



bash scripts/run\_slim\_batch.sh



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
