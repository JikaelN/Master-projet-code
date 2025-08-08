#!/bin/bash

SEED=$1
module load slim/5.0

SCRIPT_DIR=/work/FAC/FBM/DEE/jgoudet/pop_fst/jikael/slim_scripts
OUTPUT_DIR=/scratch/jntoko/trees_files/stepping
mkdir -p "$OUTPUT_DIR"

slim -d SEED=$SEED "$SCRIPT_DIR/stepping_neutral.slim"

TREE_FILE="stepping_${SEED}.trees"
if [ -f "$TREE_FILE" ]; then
    mv "$TREE_FILE" "$OUTPUT_DIR/"
fi
