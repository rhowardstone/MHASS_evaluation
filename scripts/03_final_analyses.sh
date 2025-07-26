#!/bin/bash
# MHASS Evaluation Pipeline - Final analysis with memory profiling
# Note the hardcoded subread accuracy values (random perturbation may vary these slightly with repetition)

# Create a timestamp for unique log file names
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGDIR="memory_profiling_logs"
mkdir -p "$LOGDIR"

# Define the time command format for detailed memory profiling
# Using verbose format plus custom format for comprehensive data
TIME_CMD="/usr/bin/time -v"
TIME_FORMAT="Command: %C\nElapsed Time: %E\nUser Time: %U seconds\nSystem Time: %S seconds\nCPU Usage: %P\nMax RSS (KB): %M\nAvg Total Memory (KB): %K\nAvg Data (KB): %D\nAvg Stack (KB): %p\nAvg Shared Text (KB): %X\nPage Size (bytes): %Z\nMajor Page Faults: %F\nMinor Page Faults: %R\nSwaps: %W\nFile System Inputs: %I\nFile System Outputs: %O\nVoluntary Context Switches: %w\nInvoluntary Context Switches: %c\nSignals: %k\nExit Status: %x"

echo "Memory profiling logs will be saved to: $LOGDIR"
echo ""

echo 'Three datasets:'
ls -sh /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/real/m54215_200618_132850.Q20.fastq /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/All_reads.fastq /data/shoreline/Simulator_datasets/Phylotag/real/Mock1-5.fastq
echo ''

#########  Phylotag:  #########
echo ''
echo '#-#-#-#-#-#-#-#-#   BEGINNING PHYLOTAG   #-#-#-#-#-#-#-#-#'
echo "Starting Phylotag analysis with memory profiling..."
$TIME_CMD -o "$LOGDIR/phylotag_memory_${TIMESTAMP}.log" \
    bash ./scripts/Get_ground_truth.sh \
    --reads /data/shoreline/Simulator_datasets/Phylotag/real/Mock1-5.fastq \
    --genomes_dir /data/shoreline/Simulator_datasets/Phylotag/genomes/ \
    --genome_labels /data/shoreline/Simulator_datasets/Phylotag/genome_labels.txt \
    --primers /data/shoreline/Simulator_datasets/Phylotag/phylotag_primers.txt \
    --output_dir /data/shoreline/Simulator_datasets/Phylotag/Analysis_final \
    --threads 180 \
    --abundances /data/shoreline/Simulator_datasets/Phylotag/genome_abunds.tsv \
    --barcodes /data/shoreline/Simulator_datasets/Phylotag/Dummy_barcodes_5.txt \
    --samples 5 --depth 22742 --min_length 1400 --max_length 1700 \
    --dist_threshold 50 --lognormal_mu 3.88 --lognormal_sigma 1.22 --subread_accuracy 0.6

echo "Phylotag memory profile saved to: $LOGDIR/phylotag_memory_${TIMESTAMP}.log"

#########  Zymo:  #########
echo ''
echo '#-#-#-#-#-#-#-#-#   BEGINNING ZYMO   #-#-#-#-#-#-#-#-#'
echo "Starting Zymo analysis with memory profiling..."
$TIME_CMD -o "$LOGDIR/zymo_memory_${TIMESTAMP}.log" \
    bash ./scripts/Get_ground_truth.sh \
    --dist_threshold 50 \
    --reads /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/real/m54215_200618_132850.Q20.fastq \
    --genomes_dir /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/sim/genomes/ \
    --genome_labels /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/sim/genome_labels.txt \
    --primers /data/shoreline/Simulator_datasets/Titan_primers.txt \
    --output_dir /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/Analysis_final \
    --threads 180 \
    --abundances /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/sim/amplicons/genome_abunds.tsv \
    --barcodes /data/shoreline/Simulator_datasets/Zymo_Titan_D6300/Strainid_or_nanoid_barcodes.txt \
    --samples 96 --depth 1939 --min_length 2000 --max_length 3000 \
    --dist_threshold 50 --lognormal_mu 3.08 --lognormal_sigma 0.94 --subread_accuracy 0.6

echo "Zymo memory profile saved to: $LOGDIR/zymo_memory_${TIMESTAMP}.log"

#########  ATCC:  #########
echo ''
echo '#-#-#-#-#-#-#-#-#   BEGINNING ATCC   #-#-#-#-#-#-#-#-#'
echo "Starting ATCC analysis with memory profiling..."
$TIME_CMD -o "$LOGDIR/atcc_memory_${TIMESTAMP}.log" \
    bash ./scripts/Get_ground_truth.sh \
    --reads /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/All_reads.fastq \
    --genomes_dir /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/genomes/assemblies/ \
    --genome_labels /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/genomes/genome_labels.tsv \
    --primers /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/primers.txt \
    --output_dir /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/Analysis_final \
    --threads 180 \
    --abundances /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/genomes/genome_abunds.tsv \
    --barcodes /data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/barcode_combinations.tsv \
    --samples 192 --depth 19133 --min_length 1400 --max_length 1700 \
    --dist_threshold 50 --lognormal_mu 3.88 --lognormal_sigma 1.22 --subread_accuracy 0.7

echo "ATCC memory profile saved to: $LOGDIR/atcc_memory_${TIMESTAMP}.log"

echo ""
echo "All memory profiling complete!"
echo "Memory logs saved in: $LOGDIR/"
echo ""
echo "To view a summary of all runs:"
echo "grep -E 'Maximum resident|Elapsed time|CPU' $LOGDIR/*_memory_${TIMESTAMP}.log"

# Optional: Create a combined summary file
SUMMARY_FILE="$LOGDIR/summary_${TIMESTAMP}.txt"
echo "Creating summary file: $SUMMARY_FILE"
echo "MHASS Memory Profiling Summary - $TIMESTAMP" > "$SUMMARY_FILE"
echo "=========================================" >> "$SUMMARY_FILE"
echo "" >> "$SUMMARY_FILE"

for logfile in "$LOGDIR"/*_memory_${TIMESTAMP}.log; do
    dataset=$(basename "$logfile" | cut -d'_' -f1)
    echo "Dataset: ${dataset^^}" >> "$SUMMARY_FILE"
    echo "-------------------" >> "$SUMMARY_FILE"
    grep -E "Maximum resident|Elapsed time|Percent of CPU|Major page faults|Minor page faults" "$logfile" >> "$SUMMARY_FILE"
    echo "" >> "$SUMMARY_FILE"
done

echo "Summary saved to: $SUMMARY_FILE"
