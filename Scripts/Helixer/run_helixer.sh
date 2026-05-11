#!/bin/bash

# Load the BUSCO module
module load helixer

Helixer.py --lineage invertebrate --fasta-path A_rustica_clean/GCA_964213965.2_icArhRust1.hap1.3_genomic.fna --species Arhopalus_rusticus --gff-output-path A_rustica_clean/Arhopalus_rusticus_hap1.3_GCA_964213945.2.gff3


