#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --account=jgoudet_pop_fst
#SBATCH --job-name=qst_fst_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=00:20:00
#SBATCH --output=/scratch/jntoko/logs/stepping/%x_%j.out
#SBATCH --error=/scratch/jntoko/logs/stepping/%x_%j.err

# Charger le module R
module load r-light/4.4.1

# Aller dans le dossier des scripts
cd /work/FAC/FBM/DEE/jgoudet/pop_fst/jikael/scripts

# Exécuter le script pour le fichier 1 avec 1000 bootstrap
Rscript --vanilla GW_Qst_Fst.r 1
