#!/bin/bash

#################### SLURM OPTIONS ####################
#SBATCH --partition cpu
#SBATCH --account jgoudet_pop_fst
#SBATCH --job-name array_qst_fst
#SBATCH --array=1-500%10
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 20G
#SBATCH --time 00:7:00
#SBATCH --output /scratch/jntoko/logs/%x_%j.out
#SBATCH --error /scratch/jntoko/logs/%x_%j.err

#################### LOAD MODULES ####################
module load r-light/4.4.1

#################### SETUP ####################
cd /work/FAC/FBM/DEE/jgoudet/pop_fst/jikael/scripts

Rscript --vanilla GW_Qst_Fst.r ${SLURM_ARRAY_TASK_ID}
