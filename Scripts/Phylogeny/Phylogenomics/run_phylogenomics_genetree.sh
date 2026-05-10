#!/bin/bash

# Define input and output directories
OUTPUT_DIR=./phylogenomics_beetle_tree_genetree

# Remove existing output directory if it exists
if [ -d "$OUTPUT_DIR" ]; then
    echo "Removing existing output directory..."
    rm -rf "$OUTPUT_DIR"
fi

BUSCO_phylogenomics.py -i ../busco_beetle_result -o "$OUTPUT_DIR" -t 48 --gene_trees_only
