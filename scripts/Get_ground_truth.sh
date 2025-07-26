#!/bin/bash
# MHASS Evaluation Pipeline

# Default values
PATCH_PRIMERS=false
DIST_THRESH=50
MAX_SEQS=120
LOGNORMAL_MU=3.88
LOGNORMAL_SIGMA=1.22
SUBREAD_ACCURACY=0.7
USE_CONSENSUS=false  # Default: use reference as ground truth for both
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
if [[ "$SCRIPT_DIR" == "." ]]; then
  SCRIPT_DIR=$(pwd)
fi

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --reads)
      INPUT_READS="$2"
      shift 2
      ;;
    --genomes_dir)
      INPUT_GENOMES="$2"
      shift 2
      ;;
    --genome_labels)
      GENOME_LABELS="$2"
      shift 2
      ;;
    --primers)
      INPUT_PRIMERS="$2"
      shift 2
      ;;
    --output_dir)
      OUTPUT_DIR="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
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
    --dist_threshold)
      DIST_THRESH="$2"
      shift 2
      ;;
    --max_seqs)
      MAX_SEQS="$2"
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
    --subread_accuracy)
      SUBREAD_ACCURACY="$2"
      shift 2
      ;;
    --patch-primers)
      PATCH_PRIMERS=true
      shift
      ;;
    --use-consensus)
      USE_CONSENSUS=true
      shift
      ;;
    --help)
      echo "Usage: $0 [options]"
      echo "Options:"
      echo "  --reads         Input FASTQ file with reads"
      echo "  --genomes_dir   Directory containing genome FASTA files"
      echo "  --genome_labels TSV file mapping genome filenames to IDs"
      echo "  --primers       Primer file for AmpliconHunter"
      echo "  --output_dir    Output directory for results"
      echo "  --threads       Number of threads to use"
      echo "  --abundances    TSV file with genome abundances"
      echo "  --barcodes      TSV file with barcodes for MHASS"
      echo "  --samples       Number of samples for simulation"
      echo "  --depth         Sequencing depth for simulation"
      echo "  --min_length    Minimum amplicon length"
      echo "  --max_length    Maximum amplicon length"
      echo "  --dist_threshold Maximum edit distance for consensus generation - ignored otherwise"
      echo "  --max_seqs      Maximum sequences per alignment group"
      echo "  --lognormal_mu  Lognormal distribution mu parameter (default: 3.88)"
      echo "  --lognormal_sigma Lognormal distribution sigma parameter (default: 1.22)"
      echo "  --subread_accuracy Mean subread accuracy for PBSIM simulation (default: 0.85)"
      echo "  --patch-reads   Input reads are 5'->3' and do NOT contain primers: add them"
      echo "  --use-consensus Use consensus approach for real reads (default: use reference for both)"
      echo "  --help          Display this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Set default values for optional parameters
THREADS=${THREADS:-16}
LMIN=${LMIN:-1000}
LMAX=${LMAX:-3000}

# Validate required parameters
if [[ -z "$INPUT_READS" || -z "$INPUT_GENOMES" || -z "$INPUT_PRIMERS" || -z "$OUTPUT_DIR" ]]; then
  echo "Error: Missing required parameters."
  echo "Required: --reads, --genomes_dir, --primers, --output_dir"
  exit 1
fi

# Create output directories
mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/temp"

# Add this function to the script
extract_genome_ids() {
    local genome_labels="$1"
    local amplicon_labels="$2"
    local output_file="$3"
    
    echo "Creating genome dictionary file..."
    
    # Extract unique genome IDs from labels files
    if [[ -f "$genome_labels" ]]; then
        # If genome_labels exists, extract second column
        cut -f2 "$genome_labels" | sort | uniq > "$output_file"
    elif [[ -f "$amplicon_labels" ]]; then
        # Otherwise try to extract from amplicon labels
        cut -f2 "$amplicon_labels" | sort | uniq > "$output_file"
    else
        # Create empty file if neither exists
        touch "$output_file"
        echo "Warning: Could not find genome labels, created empty dictionary"
    fi
    
    echo "Created genome dictionary with $(wc -l < "$output_file") entries"
}


# Function to generate amplicon labels from genome filenames

generate_amplicon_labels() {
  local ah_amplicons="$1"
  local genome_labels="$2"
  local output_file="$3"
  
  echo "Generating amplicon labels from genome sources..."
  
  # Create the header for the labels file with MHASS-compatible column names
  echo -e "asvid\tgenomeid" > "$output_file"
  
  # Process amplicon headers and extract source information
  grep ">" "$ah_amplicons" | sed 's/>//' | while read -r header; do
    amplicon_id="$header"
    # Extract source from AmpliconHunter header
    source_file=$(echo "$header" | grep -o 'source=[^\.]*\.[^\.]*' | sed 's/source=//')
    
    if [[ -n "$source_file" ]]; then
      # Get just the filename from the source path
      base_source=$(basename "$source_file")
      
      # If genome_labels file is provided, use it to map source to genome ID
      if [[ -n "$genome_labels" && -f "$genome_labels" ]]; then
        # Look up the clean genome ID from the genome_labels file
        genome_id=$(grep -F "$base_source" "$genome_labels" | cut -f2)
        
        if [[ -z "$genome_id" ]]; then
          # If not found, use the source filename as is
          genome_id="$base_source"
        fi
      else
        # Use the source filename as the genome ID
        genome_id="$base_source"
      fi
      
      # Add to mapping file
      echo -e "$amplicon_id\t$genome_id" >> "$output_file"
    fi
  done
  
  echo "Created amplicon labels file: $output_file"
}

echo "=== MHASS Evaluation Pipeline ==="
echo "Input reads:      $INPUT_READS"
echo "Genomes dir:      $INPUT_GENOMES"
echo "Genome labels:    $GENOME_LABELS"
echo "Primers:          $INPUT_PRIMERS"
echo "Output dir:       $OUTPUT_DIR"
echo "Threads:          $THREADS"
echo "Length range:     $LMIN-$LMAX"
echo "Distance thresh:  $DIST_THRESH"
echo "Max seqs:         $MAX_SEQS"
echo "Lognormal mu:     $LOGNORMAL_MU"
echo "Lognormal sigma:  $LOGNORMAL_SIGMA"
echo "Subread accuracy: $SUBREAD_ACCURACY"  # Add this line
echo "Use consensus:    $USE_CONSENSUS"
echo ""

echo "1) Extracting Reference Sequences"
# Find all genome files in the genomes directory
find "$INPUT_GENOMES" -type f \( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \) > "$OUTPUT_DIR/genome_filenames.txt"
echo "Found $(wc -l < "$OUTPUT_DIR/genome_filenames.txt") genome files."

# Run AmpliconHunter to extract amplicons
echo "Running AmpliconHunter to extract amplicons..."
ampliconhunter run --Lmin "$LMIN" --Lmax "$LMAX" --Tm 50 --clamp 1 --mismatches 2 \
  "$OUTPUT_DIR/genome_filenames.txt" "$INPUT_PRIMERS" "$OUTPUT_DIR/ref/"

# Filter out off-target amplicons
echo "Filtering out off-target amplicons..."
python "$SCRIPT_DIR/Filter_out_offtarget.py" "$OUTPUT_DIR/ref/amplicons.fa" "$OUTPUT_DIR/ref/ref-filtered.fa"
INPUT_REF="$OUTPUT_DIR/ref/ref-filtered.fa"

# Generate amplicon labels
LABELS="$OUTPUT_DIR/amplicon_labels.txt"
generate_amplicon_labels "$INPUT_REF" "$GENOME_LABELS" "$LABELS"

# Create a genome dictionary file with all unique genome IDs
GENOME_DICT="$OUTPUT_DIR/genome_dictionary.txt"
extract_genome_ids "$GENOME_LABELS" "$LABELS" "$GENOME_DICT"


echo "1) Dedup ref:"
python "$SCRIPT_DIR/Dedup_fasta.py" \
  -i "$INPUT_REF" \
  -o "$OUTPUT_DIR/ref-deduped.fa" \
  -m "$OUTPUT_DIR/ref-dedup-mapping.tsv" \
  -g "$GENOME_LABELS"


echo "2) Run AH on real reads:"
if [[ "$INPUT_READS" == *.fastq || "$INPUT_READS" == *.fq ]]; then
  echo "Converting FASTQ to FASTA format..."
  awk 'NR%4==1 {printf(">%s\n", substr($0,2))} NR%4==2 {print}' "$INPUT_READS" > "$OUTPUT_DIR/reads.fa"
else
  cp "$INPUT_READS" "$OUTPUT_DIR/reads.fa"
fi

if [[ "$PATCH_PRIMERS" == true ]]; then
  echo "Patching primers onto reads..."
  python "$SCRIPT_DIR/Patch_primers_onto_reads.py" \
    --reads "$OUTPUT_DIR/reads.fa" \
    --primers "$INPUT_PRIMERS" \
    --output "$OUTPUT_DIR/reads_patched.fa"
  PATCHED_READS="$OUTPUT_DIR/reads_patched.fa"
else
  PATCHED_READS="$OUTPUT_DIR/reads.fa"
fi

# Split FASTA for parallel processing
mkdir -p "$OUTPUT_DIR/reads"
python "$SCRIPT_DIR/Split_fasta.py" $PATCHED_READS -o "$OUTPUT_DIR/reads" -n $((2 * THREADS))
ls "$OUTPUT_DIR/reads/"* > "$OUTPUT_DIR/reads.txt"

# Run AmpliconHunter for read processing
echo "Running AmpliconHunter on reads..."
ampliconhunter run "$OUTPUT_DIR/reads.txt" "$INPUT_PRIMERS" "$OUTPUT_DIR/demux" \
  --mismatches 3 --clamp 1 --Tm 0 --Lmin "$LMIN" --Lmax "$LMAX"

python "$SCRIPT_DIR/Filter_out_offtarget.py" "$OUTPUT_DIR/demux/amplicons.fa" "$OUTPUT_DIR/demux/amplicons.filtered.fa"
mv "$OUTPUT_DIR/demux/amplicons.filtered.fa" "$OUTPUT_DIR/demux/amplicons.fa"


echo "3) Calculate Edit distances:"
python "$SCRIPT_DIR/Calculate_Edit_Distance.py" --input1 "$OUTPUT_DIR/ref-deduped.fa" \
  --output "$OUTPUT_DIR/demux/dists.tsv" \
  --temp_dir "$OUTPUT_DIR/temp" --input2 "$OUTPUT_DIR/demux/amplicons.fa" \
  --threads "$THREADS" --mode NW --closest_output "$OUTPUT_DIR/demux/closest.txt"

# Conditional consensus generation based on --use-consensus flag
if [[ "$USE_CONSENSUS" == true ]]; then
  echo "Using consensus approach for real reads..."
  
  echo '4) Establish ground truth for real mock community:'
  python "$SCRIPT_DIR/Group_reads_by_ground_truth_ASV.py" -c "$OUTPUT_DIR/demux/closest.txt" \
    -a "$OUTPUT_DIR/demux/amplicons.fa" \
    -o "$OUTPUT_DIR/grouped" \
    -t random -d "$DIST_THRESH"

  echo " Randomly subsample grouped sequences (max $MAX_SEQS per file):"
  mkdir -p "$OUTPUT_DIR/grouped_subsampled"
  export MAX_SEQS OUTPUT_DIR SCRIPT_DIR
  ls "$OUTPUT_DIR"/grouped/*.fa | \
    parallel -j "$THREADS" \
    'python '"$SCRIPT_DIR"'/Subsample_fasta.py \
          -i {} \
          -o '"$OUTPUT_DIR"'/grouped_subsampled/{/.}.fa \
          -n '"$MAX_SEQS"' -r'

  echo " Align with MAFFT (high quality settings), generate consensus sequences:"
  ALIGNED_DIR="$OUTPUT_DIR/grouped_aligned"
  mkdir -p "$ALIGNED_DIR"
  export ALIGNED_DIR

  # Run MAFFT in parallel
  time parallel -j "$THREADS" \
    'base=$(basename {} .fa); mafft --maxiterate 1000 --localpair --quiet {} > '"$ALIGNED_DIR"'/${base}.aligned' \
    ::: "$OUTPUT_DIR"/grouped_subsampled/*.fa

  # Generate consensus sequences
  echo ' Generate consensus sequences...'
  python "$SCRIPT_DIR/Generate_consensus_sequences.py" \
    -i "$ALIGNED_DIR" \
    -o "$OUTPUT_DIR/consensus" \
    -c "$OUTPUT_DIR/all_consensus.fa" \
    -p 0.66 -v 3 -a singular -t "$THREADS"

  echo ' Fix sparse consensus sequences (1 or 2 reads)'
  python "$SCRIPT_DIR/Fix_sparse_consensus.py" \
    -c "$OUTPUT_DIR/demux/closest.txt" \
    -a "$OUTPUT_DIR/demux/amplicons.fa" \
    -i "$OUTPUT_DIR/all_consensus.fa" \
    -o "$OUTPUT_DIR/all_consensus_fixed.fa" \
    -v

  # Move the fixed file to replace the original (lazy but hey)
  mv "$OUTPUT_DIR/all_consensus_fixed.fa" "$OUTPUT_DIR/all_consensus.fa"


  echo "5) Calculate edit distance between consensus sequences (\"ground-truth\") and AH-corrected reads:"
  python "$SCRIPT_DIR/Calculate_Edit_Distance.py" --input1 "$OUTPUT_DIR/all_consensus.fa" \
    --output "$OUTPUT_DIR/consensus-dists-real.tsv" \
    --temp_dir "$OUTPUT_DIR/temp" --input2 "$OUTPUT_DIR/demux/amplicons.fa" \
    --threads "$THREADS" --mode NW --closest_output "$OUTPUT_DIR/consensus-closest-real.txt"
  
  # Set the file to use for real data comparisons
  REAL_CLOSEST_FILE="$OUTPUT_DIR/consensus-closest-real.txt"
  REAL_GROUND_TRUTH_FASTA="$OUTPUT_DIR/all_consensus.fa"
  
else
  echo "Using reference sequences as ground truth for both real and simulated reads..."
  
  # Skip steps 4-5 (consensus generation)
  echo "4-5) Skipping consensus generation - using references directly"
  
  # The distances have already been calculated in step 3, so just copy/rename for consistency
  cp "$OUTPUT_DIR/demux/closest.txt" "$OUTPUT_DIR/ref-closest-real.txt"
  
  # Set the file to use for real data comparisons
  REAL_CLOSEST_FILE="$OUTPUT_DIR/ref-closest-real.txt"
  REAL_GROUND_TRUTH_FASTA="$OUTPUT_DIR/ref-deduped.fa"
fi


echo "6) use MHASS to simulate a comparable dataset:"
mkdir -p "$OUTPUT_DIR/simulated"
time mhass --amplicon-fasta "$INPUT_REF" --amplicon-genome-labels "$LABELS" \
--output-dir "$OUTPUT_DIR/simulated" \
--num-samples "$SAMPLES" --num-reads "$DEPTH" \
--genome-distribution "empirical:$ABUNDS" --barcode-file "$BARCODES" \
--subread-accuracy "$SUBREAD_ACCURACY" --np-distribution-type lognormal --lognormal-mu "$LOGNORMAL_MU" \
--lognormal-sigma "$LOGNORMAL_SIGMA" --np-min 2 --np-max 59 --threads "$THREADS"


echo "7) AmpliconHunter to correct the reads:"
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

python "$SCRIPT_DIR/Filter_out_offtarget.py" "$OUTPUT_DIR/simulated/demux/amplicons.fa" "$OUTPUT_DIR/simulated/demux/amplicons.filtered.fa"
mv "$OUTPUT_DIR/simulated/demux/amplicons.filtered.fa" "$OUTPUT_DIR/simulated/demux/amplicons.fa"


echo "8) Calculate edit distance between input reference sequences and simulated reads:"
python "$SCRIPT_DIR/Calculate_Edit_Distance.py" --input1 "$OUTPUT_DIR/ref-deduped.fa" \
--output "$OUTPUT_DIR/simulated/ref-dists-sim.tsv" \
--temp_dir "$OUTPUT_DIR/temp" --input2 "$OUTPUT_DIR/simulated/demux/amplicons.fa" \
--threads "$THREADS" --mode NW --closest_output "$OUTPUT_DIR/simulated/ref-closest-sim.txt"


echo "9) Plot edit distance distributions:"
python "$SCRIPT_DIR/Plot_dists_by_ASV_and_genome.py" \
-r "$REAL_CLOSEST_FILE" \
-s "$OUTPUT_DIR/simulated/ref-closest-sim.txt" \
-l "$LABELS" \
-o "$OUTPUT_DIR/plots" 


echo "10) Error pattern analysis..."
echo "  Analyzing error patterns in real data:"
python "$SCRIPT_DIR/Analyze_error_patterns.py" \
  --closest_file "$REAL_CLOSEST_FILE" \
  --reads_fasta "$OUTPUT_DIR/demux/amplicons.fa" \
  --reference_fasta "$REAL_GROUND_TRUTH_FASTA" \
  --output_dir "$OUTPUT_DIR/error_analysis_real" \
  --sample_size 100000 \
  --data_type real

echo "  Analyzing error patterns in simulated data (reads vs reference):"
python "$SCRIPT_DIR/Analyze_error_patterns.py" \
  --closest_file "$OUTPUT_DIR/simulated/ref-closest-sim.txt" \
  --reads_fasta "$OUTPUT_DIR/simulated/demux/amplicons.fa" \
  --reference_fasta "$OUTPUT_DIR/ref-deduped.fa" \
  --output_dir "$OUTPUT_DIR/error_analysis_simulated" \
  --sample_size 100000 \
  --data_type simulated


echo "11) Compare sequence length distributions:"
if [[ "$USE_CONSENSUS" == true ]]; then
  python "$SCRIPT_DIR/Compare_length_distributions.py" \
  -1 "$OUTPUT_DIR/ref-deduped.fa" \
  -2 "$OUTPUT_DIR/all_consensus.fa" \
  -n1 "Reference ASVs" \
  -n2 "Consensus ASVs" \
  -o "$OUTPUT_DIR/length_comparison" \
  --bins 50 \
  --min_length "$LMIN" \
  --max_length "$LMAX"
else
  echo "Skipping length distribution comparison (no consensus sequences generated)"
fi


echo "12) Per-genome distribution plots:"
python "$SCRIPT_DIR/Plot_nearest_dists_by_genome_and_ASV.py" \
--real_closest "$REAL_CLOSEST_FILE" \
--sim_closest "$OUTPUT_DIR/simulated/ref-closest-sim.txt" \
--outdir "$OUTPUT_DIR/per_genome_plots" 


# Steps 13-14 only make sense with consensus approach
if [[ "$USE_CONSENSUS" == true ]]; then
  echo "13) Create heatmap of reference vs consensus distances:"
  python "$SCRIPT_DIR/Heatmap_ref_vs_consensus.py" \
  --ref_fasta "$OUTPUT_DIR/ref-deduped.fa" \
  --cons_fasta "$OUTPUT_DIR/all_consensus.fa" \
  --out_png "$OUTPUT_DIR/ref_vs_consensus_heatmap.png" \
  --labels "$LABELS" \
  --threads "$THREADS" --log

  echo "14) Create all-vs-all distance heatmaps within each file:"
  python "$SCRIPT_DIR/Heatmap_all_vs_all.py" \
  --ref_fasta "$OUTPUT_DIR/ref-deduped.fa" \
  --cons_fasta "$OUTPUT_DIR/all_consensus.fa" \
  --output_dir "$OUTPUT_DIR/all_vs_all_heatmaps" \
  --labels "$LABELS" \
  --threads "$THREADS"  --log
else
  echo "13-14) Skipping consensus-specific heatmaps (no consensus sequences generated)"
fi


echo "15) Plot expected vs observed abundances (ASV-level, no threshold):"
python "$SCRIPT_DIR/Plot_expected_vs_observed.py" \
  --dedup_mapping "$OUTPUT_DIR/ref-dedup-mapping.tsv" \
  --genome_abunds "$ABUNDS" \
  --real_closest "$REAL_CLOSEST_FILE" \
  --sim_closest "$OUTPUT_DIR/simulated/ref-closest-sim.txt" \
  --out_dir "$OUTPUT_DIR/abundance_plots_asv" \
  --simplify_ids


echo "16) Plot expected vs observed abundances (Genome-level, no threshold):"
python "$SCRIPT_DIR/Plot_expected_vs_observed.py" \
  --dedup_mapping "$OUTPUT_DIR/ref-dedup-mapping.tsv" \
  --genome_abunds "$ABUNDS" \
  --real_closest "$REAL_CLOSEST_FILE" \
  --sim_closest "$OUTPUT_DIR/simulated/ref-closest-sim.txt" \
  --out_dir "$OUTPUT_DIR/abundance_plots_genome" \
  --simplify_ids \
  --summarize_by_genome
  

# Step 17 only makes sense with consensus approach
if [[ "$USE_CONSENSUS" == true ]]; then
  echo "17) Visualize alignments:"
  python "$SCRIPT_DIR/Visualize_alignments.py" \
    -i "$ALIGNED_DIR" \
    -o "$OUTPUT_DIR/alignment_visualizations" \
    -f png \
    --max_seqs 100 \
    --text_output
else
  echo "17) Skipping alignment visualizations (no alignments generated)"
fi
  
  
echo "18) Check for ASV-level misassignments in simulated data:"
python "$SCRIPT_DIR/Check_ASV_misassignments.py" \
  --closest_file "$OUTPUT_DIR/simulated/ref-closest-sim.txt" \
  --output_dir "$OUTPUT_DIR/misassignment_analysis_sim" \
  --max_distance 99999 \
  --data_type simulated \
  --min_reads 1



echo "Finally) Cleanup:"
rm -rf "$OUTPUT_DIR/reads" "$OUTPUT_DIR/simulated/corrected-reads"
if [[ "$USE_CONSENSUS" == true ]]; then
  rm -rf "$OUTPUT_DIR/consensus"
fi


echo "Analysis complete! All results saved to: $OUTPUT_DIR"
echo "Ground truth approach used: $(if [[ "$USE_CONSENSUS" == true ]]; then echo "Consensus"; else echo "Reference"; fi)"
