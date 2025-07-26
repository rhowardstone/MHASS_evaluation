#!/bin/bash

# Define source directories
ATCC_DIR="/data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/Analysis_final"
PHYLOTAG_DIR="/data/shoreline/Simulator_datasets/Phylotag/Analysis_final"
ZYMO_DIR="/data/shoreline/Simulator_datasets/Zymo_Titan_D6300/Analysis_final"

# Define target directory
UNIFIED_DIR="/data/shoreline/Simulator_datasets/Unified_results"

# Create target directory if it doesn't exist
mkdir -p "$UNIFIED_DIR"

# Create per_genome_plots subdirectories for each dataset
mkdir -p "$UNIFIED_DIR/per_genome_plots_ATCC"
mkdir -p "$UNIFIED_DIR/per_genome_plots_Phylotag"
mkdir -p "$UNIFIED_DIR/per_genome_plots_Zymo"

# Function to find and copy PNG files from a directory
find_and_copy_pngs() {
    local source_dir=$1
    local target_dir=$2
    local suffix=$3
    
    if [ -d "$source_dir" ]; then
        echo "  Processing $source_dir..."
        find "$source_dir" -name "*.png" -type f | while read -r png_file; do
            filename=$(basename "$png_file")
            base_name="${filename%.png}"
            new_name="${base_name}_${suffix}.png"
            cp -v "$png_file" "$target_dir/$new_name"
        done
    else
        echo "  Directory not found: $source_dir"
    fi
}

# Process each dataset
process_dataset() {
    local source_dir=$1
    local dataset_name=$2
    
    echo "Processing $dataset_name dataset from $source_dir..."
    
    # First, let's see what's actually in the Analysis directory
    echo "  Contents of $source_dir:"
    ls -la "$source_dir" | grep -E "(plots|analysis|error|abundance)" | head -20
    
    # Process known subdirectories
    # Based on your Get_ground_truth.sh script, these are the likely directories:
    
    # Main plots directory
    if [ -d "$source_dir/plots" ]; then
        mkdir -p "$UNIFIED_DIR/plots"
        find_and_copy_pngs "$source_dir/plots" "$UNIFIED_DIR/plots" "$dataset_name"
    fi
    
    # Error analysis directories
    if [ -d "$source_dir/error_analysis_real" ]; then
        mkdir -p "$UNIFIED_DIR/error_analysis_real"
        find_and_copy_pngs "$source_dir/error_analysis_real" "$UNIFIED_DIR/error_analysis_real" "$dataset_name"
    fi
    
    if [ -d "$source_dir/error_analysis_simulated" ]; then
        mkdir -p "$UNIFIED_DIR/error_analysis_simulated"
        find_and_copy_pngs "$source_dir/error_analysis_simulated" "$UNIFIED_DIR/error_analysis_simulated" "$dataset_name"
    fi
    
    # Abundance plots
    if [ -d "$source_dir/abundance_plots_asv" ]; then
        mkdir -p "$UNIFIED_DIR/abundance_plots_asv"
        find_and_copy_pngs "$source_dir/abundance_plots_asv" "$UNIFIED_DIR/abundance_plots_asv" "$dataset_name"
    fi
    
    if [ -d "$source_dir/abundance_plots_genome" ]; then
        mkdir -p "$UNIFIED_DIR/abundance_plots_genome"
        find_and_copy_pngs "$source_dir/abundance_plots_genome" "$UNIFIED_DIR/abundance_plots_genome" "$dataset_name"
    fi
    
    # Misassignment analysis
    if [ -d "$source_dir/misassignment_analysis_sim" ]; then
        mkdir -p "$UNIFIED_DIR/misassignment_analysis_sim"
        find_and_copy_pngs "$source_dir/misassignment_analysis_sim" "$UNIFIED_DIR/misassignment_analysis_sim" "$dataset_name"
    fi
    
    # Length comparison (if exists)
    if [ -d "$source_dir/length_comparison" ]; then
        mkdir -p "$UNIFIED_DIR/length_comparison"
        find_and_copy_pngs "$source_dir/length_comparison" "$UNIFIED_DIR/length_comparison" "$dataset_name"
    fi
    
    # Heatmaps (if exist)
    if [ -d "$source_dir/all_vs_all_heatmaps" ]; then
        mkdir -p "$UNIFIED_DIR/all_vs_all_heatmaps"
        find_and_copy_pngs "$source_dir/all_vs_all_heatmaps" "$UNIFIED_DIR/all_vs_all_heatmaps" "$dataset_name"
    fi
    
    # Direct PNG files in the Analysis directory
    echo "  Looking for PNG files directly in $source_dir..."
    find "$source_dir" -maxdepth 1 -name "*.png" -type f | while read -r png_file; do
        filename=$(basename "$png_file")
        base_name="${filename%.png}"
        new_name="${base_name}_${dataset_name}.png"
        mkdir -p "$UNIFIED_DIR/root_level_plots"
        cp -v "$png_file" "$UNIFIED_DIR/root_level_plots/$new_name"
    done
    
    # Process per_genome_plots separately (keep in dataset-specific folders)
    if [ -d "$source_dir/per_genome_plots" ]; then
        echo "  Processing per_genome_plots..."
        find "$source_dir/per_genome_plots" -name "*.png" -type f | while read -r png_file; do
            filename=$(basename "$png_file")
            cp -v "$png_file" "$UNIFIED_DIR/per_genome_plots_${dataset_name}/"
        done
    fi
    
    echo ""
}

# Process each dataset
echo "Starting to process datasets..."
echo ""

# Check if the Analysis directories exist
for dir in "$ATCC_DIR" "$PHYLOTAG_DIR" "$ZYMO_DIR"; do
    if [ ! -d "$dir" ]; then
        echo "WARNING: Directory does not exist: $dir"
        # Try without _final suffix
        alt_dir="${dir%_final}"
        if [ -d "$alt_dir" ]; then
            echo "  Found alternative: $alt_dir"
        fi
    fi
done

# Process datasets (try both with and without _final suffix)
if [ -d "$ATCC_DIR" ]; then
    process_dataset "$ATCC_DIR" "ATCC"
elif [ -d "${ATCC_DIR%_final}" ]; then
    process_dataset "${ATCC_DIR%_final}" "ATCC"
fi

if [ -d "$PHYLOTAG_DIR" ]; then
    process_dataset "$PHYLOTAG_DIR" "Phylotag"
elif [ -d "${PHYLOTAG_DIR%_final}" ]; then
    process_dataset "${PHYLOTAG_DIR%_final}" "Phylotag"
fi

if [ -d "$ZYMO_DIR" ]; then
    process_dataset "$ZYMO_DIR" "Zymo"
elif [ -d "${ZYMO_DIR%_final}" ]; then
    process_dataset "${ZYMO_DIR%_final}" "Zymo"
fi

# Process line search analysis files
echo "Processing line search analysis files..."
echo ""

# Create linesearch_analysis directory with subdirectories
mkdir -p "$UNIFIED_DIR/linesearch_analysis/ATCC"
mkdir -p "$UNIFIED_DIR/linesearch_analysis/Phylotag"
mkdir -p "$UNIFIED_DIR/linesearch_analysis/Zymo"

# Function to copy linesearch files
copy_linesearch_files() {
    local base_dir=$1
    local dataset_name=$2
    local target_subdir=$3
    
    echo "Processing linesearch files for $dataset_name..."
    
    # Define the linesearch directory
    local linesearch_dir="${base_dir}/linesearch_analysis"
    
    if [ -d "$linesearch_dir" ]; then
        # Copy subread_accuracy_effects.png
        if [ -f "$linesearch_dir/subread_accuracy_effects.png" ]; then
            cp -v "$linesearch_dir/subread_accuracy_effects.png" "$UNIFIED_DIR/linesearch_analysis/$target_subdir/"
        else
            echo "  Warning: subread_accuracy_effects.png not found"
        fi
        
        # Find and copy the best accuracy comparison plots
        # Look for files matching the pattern Analysis-acc*_comparison_plots.png
        find "$linesearch_dir" -name "Analysis-acc*_comparison_plots.png" -type f | sort | tail -2 | while read -r file; do
            if [ -f "$file" ]; then
                filename=$(basename "$file")
                cp -v "$file" "$UNIFIED_DIR/linesearch_analysis/$target_subdir/$filename"
            fi
        done
    else
        echo "  Warning: linesearch_analysis directory not found at $linesearch_dir"
    fi
    echo ""
}

# Copy linesearch files for each dataset
copy_linesearch_files "/data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003" "ATCC" "ATCC"
copy_linesearch_files "/data/shoreline/Simulator_datasets/Phylotag" "Phylotag" "Phylotag"
copy_linesearch_files "/data/shoreline/Simulator_datasets/Zymo_Titan_D6300" "Zymo" "Zymo"

echo "Files have been organized in $UNIFIED_DIR"
echo ""
echo "Summary of what was created:"
echo "------------------------"
find "$UNIFIED_DIR" -type d | sort | while read -r dir; do
    count=$(find "$dir" -name "*.png" 2>/dev/null | wc -l)
    if [ "$count" -gt 0 ]; then
        echo "$dir: $count PNG files"
    fi
done
