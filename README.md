# Stevens-element-analysis

This repository contains analysis scripts and summary tables for comparative genomics of beetles, with a focus on synteny, transposable elements, telomere repeats, and sex chromosome evolution.

## Repository structure

- `Data/` - curated summary tables and utility data used in analysis
- `Scripts/Phylogeny/Phylopairs/AC/` and `Scripts/Phylogeny/Phylopairs/IR/` - alignment chain data for Time_calibrated_cai_pipeline
- `Scripts/` - analysis pipelines and helper scripts organized by task
  - `data_preparation/` - genome extraction, filtering, and BUSCO summarization utilities
  - `annotation/` - genome annotation and prediction pipeline wrappers
    - `annotation/Busco/` - BUSCO batch execution scripts
    - `annotation/EarlGrey/` - earlGrey genome analysis wrappers
    - `annotation/Helixer/` - Helixer prediction pipeline scripts
  - `Phylogeny/` - phylogeny and tree inference workflows
    - `Phylogeny/Astral/`
    - `Phylogeny/IQTree3/`
    - `Phylogeny/Phylogenomics/`
  - `analysis/` - downstream analyses and plotting workflows
    - `analysis/Genespace/` - gene space and chromosome plotting scripts
    - `analysis/TE/` - transposable element analysis scripts
    - `analysis/Tidk/` - tidk pipelines and plotting
    - `analysis/Ancestral_X_analysis/` - X chromosome focused analyses and plots
    - `analysis/Phylopairs/`
- `LICENSE` - license for the repository

## Info

- pipeline scripts for data extraction, BUSCO summarization, phylogenetics, and plotting
- summary tables in `Data/` used for downstream figures and analysis
- example shell scripts for running jobs on a cluster or local environment

## What is not included
Some large input data and runtime outputs are not part of this repository:

- raw or downloaded genome archives, all genome files information can be found in Supplemental table 1
- extracted `.fna` genome files and genome sequence inputs
- BUSCO run outputs
- earlGrey run output directories

## Figshare repository
Figshare: https://figshare.com/projects/Datasets_for_Species-rich_and_genomically_diverse_comparative_genomics_reveal_how_fusions_fissions_and_sex_chromosomes_have_shaped_beetle_evolution_/274777
Will have Helixer annotations, BUSCO output: full_table, Earl Grey annotations, and Genespace Bed, Peptide, and orthofinder directories.

## Recommended workflow

# Phylogenetic tree
1. Prepare input genomes and place them in a desired data directory.
2. Run `Scripts/data_preparation/extract_fna.sh` to extract FASTA files from genome archives.
3. Run BUSCO on the extracted genomes with `Scripts/annotation/Busco/busco_b_parallel.sh`.
4. Run the phylogenomics pipeline from `Scripts/Phylogeny/Phylogenomics/run_phylogenomics_genetree.sh`.
5. Run ASTRAL `Scripts/Phylogeny/Astral/run_Astral.sh` on the ALL.tree files from run_phylogenomics_genetree output.
6. Run `Scripts/Phylogeny/IQTree3/iqtree3_pipeline.sh`.

# Element_vis
1. Element_vis is a first pass visual script, you just need to have the full_table.tsv generated from running BUSCO
2. Run `Scripts/analysis/vis_ALG/run_element_vis_neosex.sh` to generate bar plots of neosex genomes, you would need to manually alter the species list within the script if you want to do a different group. 

# Genespace
1. Run `Scripts/data_preparation/filter_all.py` to remove any scaffolds and Y chromsome.
2. Run `Scripts/annotation/Helixer/run_helixer_array.sh`.
3. Divide your genomes into superfamily level, and run helper scripts in `Scripts/analysis/Genespace/Utilities` to generate BED and Peptide files from your helixer annotations.
4. Run `Scripts/analysis/Genespace/run_GENESPACE.R` per superfamily
5. Run `Scripts/analysis/Genespace/extrac_chr_order.R` and `Scripts/Genespace/extract_chr_map_from_fna.R` to create csv files to replace the chromosome labels for the next step.
6. run `Scripts/analysis/Genespace/genespace_plots.R` to create custom order, and custom chromosome label plots.

# Phylopairs
1. Run `Scripts/analysis/Phylopairs/parse_blk_genomes_v4.R`, this requires you to have generated genespace repositories and give it the paths to either Tribolium castaneum or the reference species *phasedBlks.csv in the riparian directory
2. Run `Scripts/analysis/Phylopairs/Time_calibrated_cai_pipeline.R`, this requires you to have AC and IR folders which are obtained from Cai's 2022 paper, Tree_Branch_lengths.txt, Beetle_genome_families.txt. It produces a ultrametric time calibrated phylogeny, and a visual with nodes that represent divergent times in the beetle radiation.
3. Run `Scripts/analysis/Phylopairs/launch_phylopair_jobs.sh`, it is a wrapper function for `run_phylopairs_superfamily_updated_beta_final.R` where it will launch 8-10 jobs in parallel to process each of the 8 superfamily, and the Chrysomeloidea families. 

# Ancestral_X_analysis
1. Requires GENESPACE directories
2. Run `Scripts/analysis/Ancestral_X_analysis/combined_summarize_and_extract_x.py` to extract BED and peptide files associated to the X and neo-X chromosome. And store the results into Ancestral X and Neo-X respective directories
3. Run GENESPACE on the two new directories
4. Run any other analysis scripts, both boxplots only need the BED files, while the `scan_conserved_thresholds.R` requires orthofinder from GENESPACE to be ran. You would then need to point the script to Genespace_run/orthofinder/results/Orthogroups/Orthogroups.GeneCount.tsv.

# Tandem Elements (TE) repeats
1. Run `Scripts/annotation/run_earlgrey.sh` with earlgrey_genomes.txt where you give it a path to your genomic fna files.
2. Run `Scripts/analysis/TE/parse_earlgrey_highlevelcount_v4.py` to extract all highLevelCount.txt in a directory.
3. Run `Scripts/analysis/TE/calc_x_vs_autosome_repeat_categories_2.py` where it takes a directory with Earl Grey results and a genespace directory that has only X chromosome bed files.

# tidK
1. Run `Scripts/analysis/Tidk/run_tidk_array.sh` to do the initial analysis, it needs a list of fna files and their paths
2. Run `Scripts/analysis/Tidk/tidk_batch_plots.R` to generate visuals to more easily identify chromosomes that have telomeric repeats on chromosomes and their associated sequences.

## Dependencies

- bash-compatible shell
- Python 3
- R and common R packages used by plotting scripts
- BUSCO
- earlGrey
- tidK
- GENESPACE
- IQ-TREE and ASTRAL
- SLURM or another HPC scheduler for batch execution
- conda or module support for environment management


