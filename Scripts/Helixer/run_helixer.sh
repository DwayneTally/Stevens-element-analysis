#!/bin/bash
#SBATCH --job-name=Helixer
#SBATCH --output=logs/Helixer%j.out
#SBATCH --error=logs/Helixer%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dwtally@iu.edu
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --partition=gpu
#SBATCH --mem=128G
#SBATCH -A r00259

# Load the BUSCO module
module load helixer

Helixer.py --lineage invertebrate --fasta-path A_rustica_clean/GCA_964213965.2_icArhRust1.hap1.3_genomic.fna --species Arhopalus_rusticus --gff-output-path A_rustica_clean/Arhopalus_rusticus_hap1.3_GCA_964213945.2.gff3


