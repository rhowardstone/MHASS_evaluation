#!/bin/bash
# Line search for MHASS subread accuracy parameter tuning
# This script runs only steps 6-9 of the pipeline with different subread accuracy values

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --base_dir)
      BASE_DIR="$2"
      shift 2
      ;;
    --script_dir)
      SCRIPT_DIR="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --input_reads)
      INPUT_READS="$2"
      shift 2
      ;;
    --input_genomes)
      INPUT_GENOMES="$2"
      shift 2
      ;;
    --genome_labels)
      GENOME_LABELS="$2"
      shift 2
      ;;
    --input_primers)
      INPUT_PRIMERS="$2"
      shift 2
      ;;
    --abundances)
      ABUNDS="$2"
      shift 2
      ;;
    --barcodes)
      BARCODES="$2"
      shift 2
      ;;
    --samples)
      SAMPLES="$2"
      shift 2
      ;;
    --depth)
      DEPTH="$2"
      shift 2
      ;;
    --min_length)
      LMIN="$2"
      shift 2
      ;;
    --max_length)
      LMAX="$2"
      shift 2
      ;;
    --lognormal_mu)
      LOGNORMAL_MU="$2"
      shift 2
      ;;
    --lognormal_sigma)
      LOGNORMAL_SIGMA="$2"
      shift 2
      ;;
    --reference_analysis)
      REFERENCE_ANALYSIS="$2"
      shift 2
      ;;
    --accuracies)
      # Read space-separated list of accuracies
      SUBREAD_ACCURACIES=($2)
      shift 2
      ;;
    --help)
      echo "Usage: $0 [options]"
      echo "Options:"
      echo "  --base_dir           Base directory for dataset"
      echo "  --script_dir         Directory containing scripts"
      echo "  --threads            Number of threads to use"
      echo "  --input_reads        Input FASTQ file with reads"
      echo "  --input_genomes      Directory containing genome FASTA files"
      echo "  --genome_labels      TSV file mapping genome filenames to IDs"
      echo "  --input_primers      Primer file for AmpliconHunter"
      echo "  --abundances         TSV file with genome abundances"
      echo "  --barcodes           TSV file with barcodes for MHASS"
      echo "  --samples            Number of samples for simulation"
      echo "  --depth              Sequencing depth for simulation"
      echo "  --min_length         Minimum amplicon length"
      echo "  --max_length         Maximum amplicon length"
      echo "  --lognormal_mu       Lognormal distribution mu parameter"
      echo "  --lognormal_sigma    Lognormal distribution sigma parameter"
      echo "  --reference_analysis Path to reference analysis directory (steps 1-5)"
      echo "  --accuracies         Space-separated list of accuracy values to test"
      echo "  --help               Display this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Validate required parameters
required_params=(
  "BASE_DIR" "SCRIPT_DIR" "THREADS" "INPUT_READS" "INPUT_GENOMES" 
  "GENOME_LABELS" "INPUT_PRIMERS" "ABUNDS" "BARCODES" "SAMPLES" 
  "DEPTH" "LMIN" "LMAX" "LOGNORMAL_MU" "LOGNORMAL_SIGMA" 
  "REFERENCE_ANALYSIS"
)

for param in "${required_params[@]}"; do
  if [[ -z "${!param}" ]]; then
    echo "Error: Missing required parameter: --${param,,}"
    exit 1
  fi
done

echo "Creating output directory: $BASE_DIR/accuracy_optimization/"
OPTIMIZATION_DIR="$BASE_DIR/accuracy_optimization"
mkdir -p "$OPTIMIZATION_DIR"


# Set default accuracies if not provided
if [[ ${#SUBREAD_ACCURACIES[@]} -eq 0 ]]; then
  SUBREAD_ACCURACIES=(1.00 0.90 0.80 0.70 0.60 0.50)
fi

# Reference files from the original analysis (steps 1-5)
INPUT_REF="$REFERENCE_ANALYSIS/ref/ref-filtered.fa"
REF_DEDUPED="$REFERENCE_ANALYSIS/ref-deduped.fa"
CONSENSUS_CLOSEST_REAL="$REFERENCE_ANALYSIS/ref-closest-real.txt"
LABELS="$REFERENCE_ANALYSIS/amplicon_labels.txt"

# Create log file
LOG_FILE="$OPTIMIZATION_DIR/linesearch_log.txt"
echo "Starting line search at $(date)" | tee "$LOG_FILE"
echo "Dataset: $BASE_DIR" | tee -a "$LOG_FILE"
echo "Total accuracy values to test: ${#SUBREAD_ACCURACIES[@]}" | tee -a "$LOG_FILE"
echo "Fixed parameters:" | tee -a "$LOG_FILE"
echo "  Lognormal mu: $LOGNORMAL_MU" | tee -a "$LOG_FILE"
echo "  Lognormal sigma: $LOGNORMAL_SIGMA" | tee -a "$LOG_FILE"
echo "" | tee -a "$LOG_FILE"

# Counter for progress
total_runs=${#SUBREAD_ACCURACIES[@]}
current_run=0

# Loop through all subread accuracy values
for accuracy in "${SUBREAD_ACCURACIES[@]}"; do
    current_run=$((current_run + 1))
    
    # Create output directory name
    OUTPUT_DIR="$OPTIMIZATION_DIR/Analysis-acc${accuracy}_mu${LOGNORMAL_MU}_sigma${LOGNORMAL_SIGMA}"
    
    echo "[$current_run/$total_runs] Running with subread accuracy=$accuracy" | tee -a "$LOG_FILE"
    echo "Output directory: $OUTPUT_DIR" | tee -a "$LOG_FILE"
    
    # Skip if already exists and has results
    if [ -f "$OUTPUT_DIR/simulated/ref-closest-sim.txt" ] && [ -f "$OUTPUT_DIR/plots/overall_distribution.png" ]; then
        echo "  Skipping - results already exist" | tee -a "$LOG_FILE"
        continue
    fi
    
    # Create necessary directories
    mkdir -p "$OUTPUT_DIR/temp"
    
    # Copy necessary files from reference analysis (steps 1-5)
    echo "  Copying reference files..." | tee -a "$LOG_FILE"
    cp -r "$REFERENCE_ANALYSIS/ref" "$OUTPUT_DIR/"
    cp "$REFERENCE_ANALYSIS/ref-deduped.fa" "$OUTPUT_DIR/"
    cp "$REFERENCE_ANALYSIS/ref-dedup-mapping.tsv" "$OUTPUT_DIR/"
    cp "$REFERENCE_ANALYSIS/amplicon_labels.txt" "$OUTPUT_DIR/"
    
    # Step 6: Run MHASS simulation with current parameters
    echo "  Running MHASS simulation..." | tee -a "$LOG_FILE"
    mkdir -p "$OUTPUT_DIR/simulated"
    
    start_time=$(date +%s)
    
    mhass --amplicon-fasta "$INPUT_REF" --amplicon-genome-labels "$LABELS" \
        --output-dir "$OUTPUT_DIR/simulated" \
        --num-samples "$SAMPLES" --num-reads "$DEPTH" \
        --genome-distribution "empirical:$ABUNDS" --barcode-file "$BARCODES" \
        --subread-accuracy "$accuracy" --np-distribution-type lognormal \
        --lognormal-mu "$LOGNORMAL_MU" --lognormal-sigma "$LOGNORMAL_SIGMA" \
        --np-min 2 --np-max 59 --threads "$THREADS" 2>&1 | tee -a "$LOG_FILE"
    
    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo "  ERROR: MHASS failed!" | tee -a "$LOG_FILE"
        continue
    fi
    
    # Step 7: AmpliconHunter to correct the reads
    echo "  Running AmpliconHunter..." | tee -a "$LOG_FILE"
    
    # Convert to FASTA
    awk 'NR%4==1 {printf(">%s\n", substr($0,2))} NR%4==2 {print}' \
        "$OUTPUT_DIR/simulated/combined_reads.fastq" > "$OUTPUT_DIR/simulated/corrected-reads.fa"
    
    # Split FASTA
    mkdir -p "$OUTPUT_DIR/simulated/corrected-reads"
    python "$SCRIPT_DIR/Split_fasta.py" "$OUTPUT_DIR/simulated/corrected-reads.fa" \
        -o "$OUTPUT_DIR/simulated/corrected-reads" -n $((2 * THREADS))
    
    # Create filename list
    ls "$OUTPUT_DIR/simulated/corrected-reads/"* > "$OUTPUT_DIR/simulated/corrected-reads.txt"
    
    # Run AmpliconHunter
    ampliconhunter run "$OUTPUT_DIR/simulated/corrected-reads.txt" "$INPUT_PRIMERS" \
        "$OUTPUT_DIR/simulated/demux" --mismatches 3 --clamp 1 --Tm 0 --Lmin "$LMIN" --Lmax "$LMAX"
    
    python "$SCRIPT_DIR/Filter_out_offtarget.py" "$OUTPUT_DIR/simulated/demux/amplicons.fa" \
        "$OUTPUT_DIR/simulated/demux/amplicons.filtered.fa"
    mv "$OUTPUT_DIR/simulated/demux/amplicons.filtered.fa" "$OUTPUT_DIR/simulated/demux/amplicons.fa"
    
    # Step 8: Calculate edit distances
    echo "  Calculating edit distances..." | tee -a "$LOG_FILE"
    python "$SCRIPT_DIR/Calculate_Edit_Distance.py" --input1 "$OUTPUT_DIR/ref-deduped.fa" \
        --output "$OUTPUT_DIR/simulated/ref-dists-sim.tsv" \
        --temp_dir "$OUTPUT_DIR/temp" --input2 "$OUTPUT_DIR/simulated/demux/amplicons.fa" \
        --threads "$THREADS" --mode NW --closest_output "$OUTPUT_DIR/simulated/ref-closest-sim.txt"
    
    # Step 9: Plot edit distance distributions
    echo "  Creating plots..." | tee -a "$LOG_FILE"
    python "$SCRIPT_DIR/Plot_dists_by_ASV_and_genome.py" \
        -r "$CONSENSUS_CLOSEST_REAL" \
        -s "$OUTPUT_DIR/simulated/ref-closest-sim.txt" \
        -l "$LABELS" \
        -o "$OUTPUT_DIR/plots"
    
    # Cleanup temporary files to save space
    echo "  Cleaning up..." | tee -a "$LOG_FILE"
    rm -rf "$OUTPUT_DIR/simulated/corrected-reads"
    rm -f "$OUTPUT_DIR/simulated/corrected-reads.fa"
    rm -f "$OUTPUT_DIR/simulated/corrected-reads.txt"
    
    end_time=$(date +%s)
    runtime=$((end_time - start_time))
    
    echo "  Completed in $runtime seconds" | tee -a "$LOG_FILE"
    echo "" | tee -a "$LOG_FILE"
done

echo "Line search completed at $(date)" | tee -a "$LOG_FILE"

# Create summary of all runs
echo "Creating summary of all parameter combinations..." | tee -a "$LOG_FILE"
summary_file="$OPTIMIZATION_DIR/linesearch_summary.csv"
echo "accuracy,mu,sigma,output_dir,completed" > "$summary_file"

for accuracy in "${SUBREAD_ACCURACIES[@]}"; do
    output_dir="$OPTIMIZATION_DIR/Analysis-acc${accuracy}_mu${LOGNORMAL_MU}_sigma${LOGNORMAL_SIGMA}"
    if [ -f "$output_dir/simulated/ref-closest-sim.txt" ]; then
        completed="yes"
    else
        completed="no"
    fi
    echo "$accuracy,$LOGNORMAL_MU,$LOGNORMAL_SIGMA,$output_dir,$completed" >> "$summary_file"
done

echo "Summary saved to $summary_file" | tee -a "$LOG_FILE"
