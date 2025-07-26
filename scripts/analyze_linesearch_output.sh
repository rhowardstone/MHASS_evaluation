#!/bin/bash
# run_complete_linesearch_modular.sh
# Post-simulation analysis wrapper script

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dataset_name)
      DATASET_NAME="$2"
      shift 2
      ;;
    --base_dir)
      BASE_DIR="$2"
      shift 2
      ;;
    --script_dir)
      SCRIPT_DIR="$2"
      shift 2
      ;;
    --reference_analysis)
      REFERENCE_ANALYSIS="$2"
      shift 2
      ;;
    --linesearch_output)
      LINESEARCH_OUTPUT="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;;
    --skip_linesearch)
      SKIP_LINESEARCH=true
      shift
      ;;
    --help)
      echo "Usage: $0 [options]"
      echo "Options:"
      echo "  --dataset_name       Name of the dataset (e.g., Zymo, ATCC, Phylotag)"
      echo "  --base_dir           Base directory for dataset"
      echo "  --script_dir         Directory containing scripts"
      echo "  --reference_analysis Path to reference analysis directory"
      echo "  --linesearch_output  Output directory for line search analysis"
      echo "  --threads            Number of threads for analysis (default: 120)"
      echo "  --skip_linesearch    Skip running the line search (just analyze existing results)"
      echo "  --help               Display this help message"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

# Set defaults
THREADS=${THREADS:-120}
SKIP_LINESEARCH=${SKIP_LINESEARCH:-false}

# Validate required parameters
if [[ -z "$DATASET_NAME" || -z "$BASE_DIR" || -z "$SCRIPT_DIR" || -z "$REFERENCE_ANALYSIS" || -z "$LINESEARCH_OUTPUT" ]]; then
  echo "Error: Missing required parameters."
  echo "Run with --help for usage information."
  exit 1
fi

# Get the real_closest file path
REAL_CLOSEST="$REFERENCE_ANALYSIS/ref-closest-real.txt"

if [ ! -f "$REAL_CLOSEST" ]; then
  echo "Error: Could not find real_closest file at: $REAL_CLOSEST"
  exit 1
fi

echo -e "\nStep 2: Analyzing line search results for $DATASET_NAME..."
mkdir -p "$LINESEARCH_OUTPUT"

# Make sure the Python scripts are executable
chmod +x "$SCRIPT_DIR/measure_distribution_differences.py"
chmod +x "$SCRIPT_DIR/analyze_linesearch_results.py"

# Run the analysis
python "$SCRIPT_DIR/analyze_linesearch_results.py" \
    --base_dir "$BASE_DIR/accuracy_optimization" \
    --real_closest "$REAL_CLOSEST" \
    --output_dir "$LINESEARCH_OUTPUT" \
    --measure_script "$SCRIPT_DIR/measure_distribution_differences.py" \
    --threads "$THREADS"

echo -e "\nLine search analysis complete for $DATASET_NAME!"
echo "Results saved to: $LINESEARCH_OUTPUT"
echo ""
echo "Key outputs:"
echo "  - $LINESEARCH_OUTPUT/all_linesearch_results.csv - All metrics for all accuracy values"
echo "  - $LINESEARCH_OUTPUT/linesearch_summary_report.txt - Summary report with best parameters"
echo "  - $LINESEARCH_OUTPUT/subread_accuracy_effects.png - Line plots showing metric vs accuracy"
echo "  - $LINESEARCH_OUTPUT/normalized_metrics_comparison.png - Normalized comparison of all metrics"
echo ""

# Display quick summary
if [ -f "$LINESEARCH_OUTPUT/all_linesearch_results.csv" ]; then
    echo "Quick results preview for $DATASET_NAME:"
    echo "Best 3 subread accuracy values by KL divergence:"
    head -1 "$LINESEARCH_OUTPUT/all_linesearch_results.csv"
    tail -n +2 "$LINESEARCH_OUTPUT/all_linesearch_results.csv" | sort -t',' -k4,4n | head -3
fi
