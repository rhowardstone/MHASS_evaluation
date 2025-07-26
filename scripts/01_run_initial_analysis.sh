#!/bin/bash
# MHASS Evaluation Pipeline - Initial analysis
# Note the only calculations used from this in the subread accuracy optimization workflow are the real -> reference edit distances.
#  	The remaining plots were created as exploratory analyses, and this initial run kept for syntactic ease.

echo 'Three datasets:'
ls -sh /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/real/m54215_200618_132850.Q20.fastq /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/All_reads.fastq /data/shoreline/Simulator_datasets/Phylotag/real/Mock1-5.fastq
echo ''


#########  Phylotag:  #########
echo ''
echo '#-#-#-#-#-#-#-#-#   BEGINNING PHYLOTAG   #-#-#-#-#-#-#-#-#'

time bash ./scripts/Get_ground_truth.sh   --reads /data/shoreline/Simulator_datasets/Phylotag/real/Mock1-5.fastq \
   --genomes_dir /data/shoreline/Simulator_datasets/Phylotag/genomes/ \
   --genome_labels /data/shoreline/Simulator_datasets/Phylotag/genome_labels.txt \
   --primers /data/shoreline/Simulator_datasets/Phylotag/phylotag_primers.txt \
   --output_dir /data/shoreline/Simulator_datasets/Phylotag/Analysis_init \
   --threads 180 \
   --abundances /data/shoreline/Simulator_datasets/Phylotag/genome_abunds.tsv \
   --barcodes /data/shoreline/Simulator_datasets/Phylotag/Dummy_barcodes_5.txt \
   --samples 5   --depth 22742   --min_length 1400   --max_length 1700 --dist_threshold 50  --lognormal_mu 3.88 --lognormal_sigma 1.22


#########  Zymo:  #########
echo ''
echo '#-#-#-#-#-#-#-#-#   BEGINNING ZYMO   #-#-#-#-#-#-#-#-#'

time bash ./scripts/Get_ground_truth.sh   --dist_threshold 50 --reads /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/real/m54215_200618_132850.Q20.fastq \
   --genomes_dir /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/sim/genomes/ \
   --genome_labels /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/sim/genome_labels.txt \
   --primers /data/shoreline/Simulator_datasets/Titan_primers.txt \
   --output_dir /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/Analysis_init \
   --threads 180 \
   --abundances /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/sim/amplicons/genome_abunds.tsv \
   --barcodes /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/Strainid_or_nanoid_barcodes.txt \
   --samples 88   --depth 2116   --min_length 2000 --max_length 3000 --dist_threshold 50   --lognormal_mu 3.08 --lognormal_sigma 0.94


#########  ATCC:  #########
echo ''
echo '#-#-#-#-#-#-#-#-#   BEGINNING ATCC   #-#-#-#-#-#-#-#-#'

time bash ./scripts/Get_ground_truth.sh   --reads /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/All_reads.fastq \
   --genomes_dir /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/genomes/assemblies/ \
   --genome_labels /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/genomes/genome_labels.tsv \
   --primers /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/primers.txt \
   --output_dir /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/Analysis_init \
   --threads 180 \
   --abundances /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/genomes/genome_abunds.tsv \
   --barcodes /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/barcode_combinations.tsv \
   --samples 192   --depth 19133   --min_length 1400   --max_length 1700  --dist_threshold 50   --lognormal_mu 3.88 --lognormal_sigma 1.22




