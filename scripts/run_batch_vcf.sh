#!/bin/bash

#################### SLURM OPTIONS ####################
#SBATCH --partition cpu
#SBATCH --account jgoudet_pop_fst
#SBATCH --job-name vcf_batch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 14
#SBATCH --mem 50GB
#SBATCH --time 00:45:00
#SBATCH --output /scratch/jntoko/logs/stepping/vcf_batch_%j.out
#SBATCH --error /scratch/jntoko/logs/stepping/vcf_batch_%j.err

#################### LOAD MODULES ####################
module load parallel
# Activate conda environment
#eval "$(/work/FAC/FBM/DEE/jgoudet/pop_fst/jikael/envs/slim_py_env/bin/conda shell.bash hook)"
#conda activate /work/FAC/FBM/DEE/jgoudet/pop_fst/jikael/envs/slim_py_env




#################### VARIABLES ####################
SCRIPT_DIR=/work/FAC/FBM/DEE/jgoudet/pop_fst/jikael/scripts
TREE_DIR=/scratch/jntoko/trees_files/stepping
VCF_DIR=/scratch/jntoko/vcf_output/stepping

mkdir -p "$VCF_DIR"
mkdir -p /scratch/jntoko/logs

#################### PROCESS FILES ####################
# List all .trees files and launch jobs in parallel
find "$TREE_DIR" -name "*.trees" | parallel -j 14 "$SCRIPT_DIR/run_one_vcf.sh {} $VCF_DIR"
