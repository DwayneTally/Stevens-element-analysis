#!/bin/bash
#SBATCH --job-name=ASTRAL_consensus
#SBATCH --output=./logs/ASTRAL_consensus%j.out
#SBATCH --error=./logs/ASTRAL_consensus%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dwtally@iu.edu
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=general
#SBATCH --mem=80G
#SBATCH -A r00259


java -jar astral.5.7.8.jar -i ../phylogenomics_beetle_tree_genetree/gene_trees_single_copy/ALL.tree -o consensus.tre


