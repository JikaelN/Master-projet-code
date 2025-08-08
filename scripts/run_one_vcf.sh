#!/bin/bash


# Chargement de l’environnement conda



TREE_FILE=$1
VCF_OUTPUT_DIR=/scratch/jntoko/vcf_output/stepping

mkdir -p "$VCF_OUTPUT_DIR"

#python /work/FAC/FBM/DEE/jgoudet/pop_fst/jikael/scripts/process_stepping_tree.py "$TREE_FILE" "$VCF_OUTPUT_DIR"

/work/FAC/FBM/DEE/jgoudet/pop_fst/jikael/envs/slim_py_env/bin/python \
    /work/FAC/FBM/DEE/jgoudet/pop_fst/jikael/scripts/process_stepping_tree.py \
    $TREE_FILE $VCF_OUTPUT_DIR
