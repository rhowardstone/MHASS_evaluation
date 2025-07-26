#!/bin/bash
# run_all_linesearches.sh
# Actual main driver script for subread accuracy optimization
# Runs line searches for all three datasets: Zymo, ATCC, and Phylotag

# Common parameters
SCRIPT_DIR="/data/shoreline/Simulator_datasets/scripts"
THREADS=180

# Define accuracy values to test (can customize per dataset if needed)
DEFAULT_ACCURACIES="1.00 0.90 0.80 0.70 0.60 0.50"  #"0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00"

echo "=== Running line searches for all three datasets: Phylotag, then Zymo, then ATCC ==="
echo ""




#########  Phylotag Dataset  #########
echo '#-#-#-#-#-#-#-#-#   PHYLOTAG subread accuracy optimization   #-#-#-#-#-#-#-#-#'
PHYLOTAG_BASE="/data/shoreline/Simulator_datasets/Phylotag"



bash "$SCRIPT_DIR/linesearch_mhass.sh" \
    --base_dir "$PHYLOTAG_BASE" \
    --script_dir "$SCRIPT_DIR" \
    --threads "$THREADS" \
    --input_reads "$PHYLOTAG_BASE/real/Mock1-5.fastq" \
    --input_genomes "$PHYLOTAG_BASE/genomes/" \
    --genome_labels "$PHYLOTAG_BASE/genome_labels.txt" \
    --input_primers "$PHYLOTAG_BASE/phylotag_primers.txt" \
    --abundances "$PHYLOTAG_BASE/genome_abunds.tsv" \
    --barcodes "$PHYLOTAG_BASE/Dummy_barcodes_5.txt" \
    --samples 5 \
    --depth 22742 \
    --min_length 1400 \
    --max_length 1700 \
    --lognormal_mu 3.88 \
    --lognormal_sigma 1.22 \
    --reference_analysis "$PHYLOTAG_BASE/Analysis_init/" \
    --accuracies "$DEFAULT_ACCURACIES"
   



# Analyze results
bash "$SCRIPT_DIR/analyze_linesearch_output.sh" \
    --dataset_name "Phylotag" \
    --base_dir "$PHYLOTAG_BASE" \
    --script_dir "$SCRIPT_DIR" \
    --reference_analysis "$PHYLOTAG_BASE/Analysis_init/" \
    --linesearch_output "$PHYLOTAG_BASE/linesearch_analysis" \
    --threads 120 \
    --skip_linesearch

echo ""
echo "=== Phylotag subread accuracy optimization complete ==="
echo ""




#########  Zymo Dataset  #########
echo '#-#-#-#-#-#-#-#-#   ZYMO subread accuracy optimization   #-#-#-#-#-#-#-#-#'
ZYMO_BASE="/data/shoreline/Simulator_datasets/Zymo_Titan_D6300"





bash "$SCRIPT_DIR/linesearch_mhass.sh" \
    --base_dir "$ZYMO_BASE" \
    --script_dir "$SCRIPT_DIR" \
    --threads "$THREADS" \
    --input_reads "$ZYMO_BASE/real/m54215_200618_132850.Q20.fastq" \
    --input_genomes "$ZYMO_BASE/sim/genomes/" \
    --genome_labels "$ZYMO_BASE/sim/genome_labels.txt" \
    --input_primers "/data/shoreline/Simulator_datasets/Titan_primers.txt" \
    --abundances "$ZYMO_BASE/sim/amplicons/genome_abunds.tsv" \
    --barcodes "$ZYMO_BASE/Strainid_or_nanoid_barcodes.txt" \
    --samples 88 \
    --depth 2116 \
    --min_length 2000 \
    --max_length 3000 \
    --lognormal_mu 3.08 \
    --lognormal_sigma 0.94 \
    --reference_analysis "$ZYMO_BASE/Analysis_init" \
    --accuracies "$DEFAULT_ACCURACIES"






# Analyze results
bash "$SCRIPT_DIR/analyze_linesearch_output.sh" \
    --dataset_name "Zymo" \
    --base_dir "$ZYMO_BASE" \
    --script_dir "$SCRIPT_DIR" \
    --reference_analysis "$ZYMO_BASE/Analysis_init" \
    --linesearch_output "$ZYMO_BASE/linesearch_analysis" \
    --threads 120 \
    --skip_linesearch

echo ""
echo "=== Zymo subread accuracy optimization complete ==="
echo ""





#########  ATCC Dataset  #########
echo '#-#-#-#-#-#-#-#-#   ATCC subread accuracy optimization   #-#-#-#-#-#-#-#-#'
ATCC_BASE="/data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003"





bash "$SCRIPT_DIR/linesearch_mhass.sh" \
    --base_dir "$ATCC_BASE" \
    --script_dir "$SCRIPT_DIR" \
    --threads "$THREADS" \
    --input_reads "$ATCC_BASE/All_reads.fastq" \
    --input_genomes "$ATCC_BASE/genomes/assemblies/" \
    --genome_labels "$ATCC_BASE/genomes/genome_labels.tsv" \
    --input_primers "$ATCC_BASE/primers.txt" \
    --abundances "$ATCC_BASE/genomes/genome_abunds.tsv" \
    --barcodes "$ATCC_BASE/barcode_combinations.tsv" \
    --samples 192 \
    --depth 19133 \
    --min_length 1400 \
    --max_length 1700 \
    --lognormal_mu 3.88 \
    --lognormal_sigma 1.22 \
    --reference_analysis "$ATCC_BASE/Analysis_init" \
    --accuracies "$DEFAULT_ACCURACIES"




# Analyze results
bash "$SCRIPT_DIR/analyze_linesearch_output.sh" \
    --dataset_name "ATCC" \
    --base_dir "$ATCC_BASE" \
    --script_dir "$SCRIPT_DIR" \
    --reference_analysis "$ATCC_BASE/Analysis_init" \
    --linesearch_output "$ATCC_BASE/linesearch_analysis" \
    --threads 120 \
    --skip_linesearch

echo ""
echo "=== ATCC subread accuracy optimization complete ==="
echo ""


echo "=== All line searches complete! ==="
echo ""
echo "Summary locations:"
echo "  Zymo:     $ZYMO_BASE/accuracy_optimization/"
echo "  ATCC:     $ATCC_BASE/accuracy_optimization/"
echo "  Phylotag: $PHYLOTAG_BASE/accuracy_optimization/"









exit









