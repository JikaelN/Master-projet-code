#!/bin/bash
#SBATCH --partition=cpu
#SBATCH --account=jgoudet_pop_fst
#SBATCH --job-name=slim_array
#SBATCH --array=1-500%14
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/scratch/jntoko/logs/stepping/slim_array_%A_%a.out
#SBATCH --error=/scratch/jntoko/logs/stepping/slim_array_%A_%a.err

module load slim/5.0

# Fichiers et dossiers
SCRIPT_DIR=/work/FAC/FBM/DEE/jgoudet/pop_fst/jikael/scripts
SEEDS_FILE=/work/FAC/FBM/DEE/jgoudet/pop_fst/jikael/seeds.txt
OUTPUT_DIR=/scratch/jntoko/trees_files/stepping
LOG_DIR=/scratch/jntoko/logs/stepping

# Créer les répertoires si nécessaires
mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

# Récupérer le SEED depuis le fichier (correspondant à SLURM_ARRAY_TASK_ID)
SEED=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SEEDS_FILE")

# Lancer le script avec le SEED
bash "$SCRIPT_DIR/run_one_seed.sh" "$SEED"

