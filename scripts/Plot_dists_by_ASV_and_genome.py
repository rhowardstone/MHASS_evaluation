#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import re

def parse_args():
    parser = argparse.ArgumentParser(description='Plot edit distance distributions from real and simulated data.')
    parser.add_argument('-r', '--real_closest', required=True,
                        help='Closest output file for real data')
    parser.add_argument('-s', '--sim_closest', required=True,
                        help='Closest output file for simulated data')
    parser.add_argument('-l', '--labels', required=False,
                        help='Optional file with amplicon to genome labels')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Output directory for plots')
    return parser.parse_args()

def extract_genome(seq_id):
    """Extract genome from deduped ID format like 'Escherichia_coli-1-2-3'."""
    # Find the position of the first dash followed by a digit
    match = re.search(r'^(.+?)-\d', seq_id)
    if match:
        return match.group(1)
    return seq_id

def extract_indices(seq_id):
    """Extract indices from deduped ID format like 'Escherichia_coli-1-2-3'."""
    # Extract all the numeric indices after the genome name
    match = re.search(r'-(\d+(?:-\d+)*)$', seq_id)
    if match:
        return match.group(1)
    return None

def read_closest(closest_file):
    """Read the closest.txt file to get edit distances."""
    data = []
    with open(closest_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            read_id = parts[0]
            edit_distance = int(parts[1])
            closest_refs = parts[2:]
            for ref in closest_refs:
                data.append((read_id, ref, edit_distance))
    return pd.DataFrame(data, columns=['read_id', 'reference_id', 'edit_distance'])

def read_labels(labels_file):
    """Read the genome labels file if provided."""
    if not labels_file:
        return {}
    
    labels = {}
    with open(labels_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                amplicon_id = parts[0]
                genome = parts[1]
                labels[amplicon_id] = genome
    return labels

def create_consistent_bins(real_distances, sim_distances, n_bins=50):
    """Create consistent bin edges for both datasets."""
    if len(real_distances) > 0 and len(sim_distances) > 0:
        # Use full range of both datasets
        min_val = min(real_distances.min(), sim_distances.min())
        max_val = max(real_distances.max(), sim_distances.max())
    elif len(real_distances) > 0:
        min_val = real_distances.min()
        max_val = real_distances.max()
    elif len(sim_distances) > 0:
        min_val = sim_distances.min()
        max_val = sim_distances.max()
    else:
        return np.linspace(0, 50, n_bins + 1)  # Default range
    
    # Add small padding to ensure all data is included
    range_padding = (max_val - min_val) * 0.02
    min_val -= range_padding
    max_val += range_padding
    
    return np.linspace(min_val, max_val, n_bins + 1)

def create_plots(real_df, sim_df, labels, output_dir):
    """Create various plots comparing real and simulated data."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Add columns for ASV ID and genome
    def add_metadata(df):
        # The reference_id IS the ASV ID - no conversion needed
        df['asv_id'] = df['reference_id']
        
        # Extract genome from ASV ID
        df['genome'] = df['reference_id'].apply(extract_genome)
        
        # Extract indices for display
        df['indices'] = df['reference_id'].apply(extract_indices)
        
        return df
    
    real_df = add_metadata(real_df)
    sim_df = add_metadata(sim_df)
    
    # Extract distances for consistent binning
    real_distances = real_df['edit_distance'] if len(real_df) > 0 else pd.Series([], dtype=int)
    sim_distances = sim_df['edit_distance'] if len(sim_df) > 0 else pd.Series([], dtype=int)
    
    # Create consistent bins for all plots
    bins = create_consistent_bins(real_distances, sim_distances, n_bins=50)
    
    # Check if we have data from both sources
    has_real = len(real_df) > 0
    has_sim = len(sim_df) > 0
    
    print(f"\n=== BINNING INFO ===")
    print(f"Real data: {len(real_df)} reads")
    if has_real:
        print(f"  Distance range: {real_distances.min()}-{real_distances.max()}")
    print(f"Simulated data: {len(sim_df)} reads")
    if has_sim:
        print(f"  Distance range: {sim_distances.min()}-{sim_distances.max()}")
    print(f"Consistent bins: {bins[0]:.1f} to {bins[-1]:.1f} ({len(bins)-1} bins)")
    
    # 1. Overall distribution comparison - FIXED BINNING AND NORMALIZATION
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Version 1: Probability (0-1 scale)
    ax = axes[0]
    if has_real:
        weights_real = np.ones(len(real_df)) / len(real_df)
        ax.hist(real_distances, bins=bins, alpha=0.6, label=f'Real data (n={len(real_df)})', 
                color='blue', weights=weights_real)
    if has_sim:
        weights_sim = np.ones(len(sim_df)) / len(sim_df)
        ax.hist(sim_distances, bins=bins, alpha=0.6, label=f'Simulated data (n={len(sim_df)})', 
                color='red', weights=weights_sim)
    
    ax.set_xlabel("Edit Distance")
    ax.set_ylabel("Probability")
    ax.set_title("Probability Distribution")
    if has_real or has_sim:
        ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Version 2: Percentage (0-100 scale)
    ax = axes[1]
    if has_real:
        weights_real = np.ones(len(real_df)) / len(real_df) * 100
        ax.hist(real_distances, bins=bins, alpha=0.6, label=f'Real data (n={len(real_df)})', 
                color='blue', weights=weights_real)
    if has_sim:
        weights_sim = np.ones(len(sim_df)) / len(sim_df) * 100
        ax.hist(sim_distances, bins=bins, alpha=0.6, label=f'Simulated data (n={len(sim_df)})', 
                color='red', weights=weights_sim)
    
    ax.set_xlabel("Edit Distance")
    ax.set_ylabel("Percentage of Reads")
    ax.set_title("Percentage Distribution")
    if has_real or has_sim:
        ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Version 3: Raw counts (for absolute comparison)
    ax = axes[2]
    if has_real:
        ax.hist(real_distances, bins=bins, alpha=0.6, label=f'Real data (n={len(real_df)})', 
                color='blue')
    if has_sim:
        ax.hist(sim_distances, bins=bins, alpha=0.6, label=f'Simulated data (n={len(sim_df)})', 
                color='red')
    
    ax.set_xlabel("Edit Distance")
    ax.set_ylabel("Count")
    ax.set_title("Raw Count Distribution")
    if has_real or has_sim:
        ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.suptitle("Edit Distance Distribution: Real vs. Simulated Data", fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "overall_distribution.png"), dpi=300)
    plt.close()
    
    # 2. Box plots by ASV
    plt.figure(figsize=(20, 10))
    
    # Combine datasets with a source column
    real_df['source'] = 'Real'
    sim_df['source'] = 'Simulated'
    combined = pd.concat([real_df, sim_df])
    
    # Get top ASVs by occurrence
    asv_counts = combined['asv_id'].value_counts()
    top_asv_ids = asv_counts.head(20).index.tolist()
    
    # Filter for top ASVs
    plot_data = combined[combined['asv_id'].isin(top_asv_ids)]
    
    if len(plot_data) > 0:
        # Check if we have data from both sources
        unique_sources = plot_data['source'].unique()
        if len(unique_sources) > 1:
            # Create box plot with hue
            sns.boxplot(x='asv_id', y='edit_distance', hue='source', data=plot_data, palette=['blue', 'red'])
            plt.legend(title='Data Source')
        else:
            # Single source, no hue needed
            source_name = unique_sources[0]
            color = 'blue' if source_name == 'Real' else 'red'
            sns.boxplot(x='asv_id', y='edit_distance', data=plot_data, color=color)
            plt.title(f"Edit Distance by ASV: {source_name} Data Only")
        
        plt.xlabel("ASV ID")
        plt.ylabel("Edit Distance")
        plt.xticks(rotation=45, ha='right')
    else:
        plt.text(0.5, 0.5, 'No data to plot', ha='center', va='center', transform=plt.gca().transAxes)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "distance_by_asv.png"), dpi=300)
    plt.close()
    
    # 3. Plot by genome with consistent binning for violin plots
    if 'genome' in combined.columns and combined['genome'].notna().any():
        plt.figure(figsize=(16, 8))
        
        # Get genomes with data
        genome_counts = combined['genome'].value_counts()
        top_genomes = genome_counts.head(10).index.tolist()
        genome_data = combined[combined['genome'].isin(top_genomes)]
        
        if len(genome_data) > 0:
            unique_sources = genome_data['source'].unique()
            if len(unique_sources) > 1:
                # Create violin plot with split
                sns.violinplot(x='genome', y='edit_distance', hue='source', data=genome_data, 
                               split=True, inner='quart', palette=['blue', 'red'])
                plt.legend(title='Data Source')
            else:
                # Single source
                source_name = unique_sources[0]
                color = 'blue' if source_name == 'Real' else 'red'
                sns.violinplot(x='genome', y='edit_distance', data=genome_data, 
                              inner='quart', color=color)
                plt.title(f"Edit Distance by Genome: {source_name} Data Only")
            
            plt.xlabel("Genome")
            plt.ylabel("Edit Distance")
            plt.xticks(rotation=45, ha='right')
        else:
            plt.text(0.5, 0.5, 'No genome data to plot', ha='center', va='center', transform=plt.gca().transAxes)
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "distance_by_genome.png"), dpi=300)
        plt.close()
        
        # 4. Heatmap of average distances by genome and ASV
        if len(real_df) > 0 and 'genome' in real_df.columns:
            plt.figure(figsize=(20, 14))
            
            try:
                # Use full ASV IDs, not just indices
                pivot_real = real_df.pivot_table(
                    index='genome', 
                    columns='asv_id',  # Use full ASV ID
                    values='edit_distance', 
                    aggfunc='mean'
                )
                
                if not pivot_real.empty and pivot_real.shape[0] > 0 and pivot_real.shape[1] > 0:
                    # Limit to reasonable number of ASVs for visibility
                    # Sort by total reads per ASV
                    asv_totals = real_df.groupby('asv_id').size().sort_values(ascending=False)
                    top_asvs = asv_totals.head(50).index  # Top 50 ASVs by read count
                    
                    # Filter pivot table to top ASVs
                    pivot_subset = pivot_real[pivot_real.columns.intersection(top_asvs)]
                    
                    # Sort columns by genome then by ASV ID
                    sorted_columns = sorted(pivot_subset.columns, key=lambda x: (extract_genome(x), x))
                    pivot_subset = pivot_subset[sorted_columns]
                    
                    # Create heatmap
                    plt.figure(figsize=(20, 12))
                    sns.heatmap(
                        pivot_subset, 
                        annot=True, 
                        fmt='.0f',  # No decimal places
                        cmap="viridis",
                        cbar_kws={'label': 'Average Edit Distance'},
                        xticklabels=True,
                        yticklabels=True
                    )
                    
                    plt.title("Average Edit Distance by Genome and ASV (Real Data)\nTop 50 ASVs by read count")
                    plt.xlabel("ASV ID")
                    plt.ylabel("Genome")
                    plt.xticks(rotation=90, ha='right')
                    plt.tight_layout()
                    plt.savefig(os.path.join(output_dir, "heatmap_genome_asv.png"), dpi=300, bbox_inches='tight')
                    
                    # Also save the data as CSV for verification
                    pivot_subset.to_csv(os.path.join(output_dir, "heatmap_genome_asv_data.csv"))
                    
                else:
                    print("Warning: Pivot table is empty, skipping heatmap")
            except Exception as e:
                print(f"Warning: Could not create heatmap: {e}")
            plt.close()
    
    # 5. Correlation between real and simulated distances
    if len(real_df) > 0 and len(sim_df) > 0:
        # Group by ASV and calculate mean distance for each dataset
        real_mean = real_df.groupby('asv_id')['edit_distance'].mean()
        sim_mean = sim_df.groupby('asv_id')['edit_distance'].mean()
        
        # Find common ASVs
        common_asvs = set(real_mean.index) & set(sim_mean.index)
        
        if common_asvs:
            # Create comparison dataframe
            comparison = pd.DataFrame({
                'real_mean_distance': real_mean[list(common_asvs)],
                'sim_mean_distance': sim_mean[list(common_asvs)]
            }).reset_index()
            
            # Create scatter plot
            plt.figure(figsize=(10, 8))
            sns.scatterplot(data=comparison, x='real_mean_distance', y='sim_mean_distance', s=80, alpha=0.7)
            
            # Add reference line
            max_val = max(comparison['real_mean_distance'].max(), comparison['sim_mean_distance'].max()) * 1.1
            plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
            
            # Add labels to points (limit to avoid overcrowding)
            for i, row in comparison.head(20).iterrows():
                plt.text(row['real_mean_distance'], row['sim_mean_distance'], 
                         row['asv_id'], fontsize=8)
            
            plt.xlabel("Mean Edit Distance (Real)")
            plt.ylabel("Mean Edit Distance (Simulated)")
            plt.title("Correlation of Edit Distances: Real vs. Simulated")
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "distance_correlation.png"), dpi=300)
            plt.close()
    
    # 6. Statistical summary
    summary_data = {
        'Real': [real_df['edit_distance'].mean() if len(real_df) > 0 else np.nan,
                 real_df['edit_distance'].median() if len(real_df) > 0 else np.nan,
                 real_df['edit_distance'].std() if len(real_df) > 0 else np.nan,
                 real_df['edit_distance'].min() if len(real_df) > 0 else np.nan,
                 real_df['edit_distance'].max() if len(real_df) > 0 else np.nan],
        'Simulated': [sim_df['edit_distance'].mean() if len(sim_df) > 0 else np.nan,
                     sim_df['edit_distance'].median() if len(sim_df) > 0 else np.nan,
                     sim_df['edit_distance'].std() if len(sim_df) > 0 else np.nan,
                     sim_df['edit_distance'].min() if len(sim_df) > 0 else np.nan,
                     sim_df['edit_distance'].max() if len(sim_df) > 0 else np.nan]
    }
    
    summary = pd.DataFrame(summary_data, index=['Mean', 'Median', 'Std Dev', 'Min', 'Max'])
    
    # Save summary to CSV
    summary.to_csv(os.path.join(output_dir, "distance_summary_stats.csv"))
    
    # Print summary
    print("\nEdit Distance Summary Statistics:")
    print(summary)
    print("\nPlots saved to:", output_dir)
    
    return summary

def main():
    args = parse_args()
    
    # Read input files
    print("Reading closest files...")
    real_df = read_closest(args.real_closest)
    sim_df = read_closest(args.sim_closest)
    
    print("Reading labels file (if provided)...")
    labels = read_labels(args.labels)
    
    # Create plots
    print("Creating plots...")
    try:
        create_plots(real_df, sim_df, labels, args.output_dir)
    except Exception as e:
        print(f"Error creating plots: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    main()
