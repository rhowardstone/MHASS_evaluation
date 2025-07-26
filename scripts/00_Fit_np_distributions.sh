#!/bin/bash

# Define paths
TITAN_PFIZER_READS=/PFIZER/30-406841736/CCS
KINNEX_READS=/data/shoreline/Simulator_datasets/NP_analysis/kinnex_16S
SCRIPTS_DIR=/data/shoreline/Simulator_datasets/scripts
OUTPUT_DIR=/data/shoreline/Simulator_datasets/NP_analysis

# Count NP distributions
echo "Extracting NP distributions..."

# Pfizer/Titan:
python $SCRIPTS_DIR/Count_np_distribution.py -j 24 --np-summary $OUTPUT_DIR/pfizer_withmax_summary.txt --minlen 2000 --maxlen 3000 --max-np 60 $TITAN_PFIZER_READS $OUTPUT_DIR/pfizer_withmax_full_summary.txt

# Kinnex/16S:
python $SCRIPTS_DIR/Count_np_distribution_bam.py -j 24 --np-summary $OUTPUT_DIR/kinnex_withmax_summary.txt --minlen 1400 --maxlen 1700 --max-np 60 $KINNEX_READS $OUTPUT_DIR/kinnex_withmax_full_summary.txt

# Analyze with new combined script (using s16 instead of 16s)
echo "Fitting distributions and creating plots..."
python $SCRIPTS_DIR/Analyze_np_distribution.py --titan_counts $OUTPUT_DIR/pfizer_withmax_summary.txt --s16_counts $OUTPUT_DIR/kinnex_withmax_summary.txt --titan_full $OUTPUT_DIR/pfizer_withmax_full_summary.txt --s16_full $OUTPUT_DIR/kinnex_withmax_full_summary.txt --out $OUTPUT_DIR/np_distribution_comparison.png

echo "Analysis complete!"
