#!/usr/bin/env python3
"""
Plot_expected_vs_observed.py
Create plots comparing expected vs observed abundances of ASVs in real and simulated data.

Expected abundance calculation:
- Each ASV (de-duped) gets a share of the genome abundance from genome_abunds.tsv
- Share is proportional to its copy number (i.e., how many identical ASVs it represents)
- If ASV 'Escherichia_coli-1' has 1 copy and 'Escherichia_coli-2/3/4' has 3 copies,
  and E. coli has 14% abundance, then expected abundances are 3.5% and 10.5% respectively.

Observed abundance:
- For real data: counts from closest file (consensus-closest-real.txt)
- For simulated data: counts from ref-closest-sim.txt
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Optional


def parse_args():
    parser = argparse.ArgumentParser(description="Create expected vs observed abundance plots.")
    parser.add_argument("--dedup_mapping", required=True, 
                        help="Deduplication mapping file (ref-dedup-mapping.tsv)")
    parser.add_argument("--genome_abunds", required=True, 
                        help="Genome abundances file (genome_abunds.tsv)")
    parser.add_argument("--real_closest", required=True, 
                        help="Closest output file for real data (consensus-closest-real.txt)")
    parser.add_argument("--sim_closest", required=True, 
                        help="Closest output file for simulated data (ref-closest-sim.txt)")
    parser.add_argument("--out_dir", required=True, 
                        help="Output directory for plots")
    parser.add_argument("--simplify_ids", action="store_true",
                        help="Convert coordinate-based IDs to simple {genome}-{index} format")
    parser.add_argument("--real_closest_threshold", type=float, default=None,
                        help="Maximum edit distance threshold for real data (optional)")
    parser.add_argument("--sim_closest_threshold", type=float, default=None,
                        help="Maximum edit distance threshold for simulated data (optional)")
    parser.add_argument("--summarize_by_genome", action="store_true",
                        help="Summarize results by genome (one point per genome)")
    
    return parser.parse_args()


def read_dedup_mapping(file_path: str) -> Dict:
    """
    Read the deduplication mapping file.
    Returns:
    - asv_to_genome: Dict mapping ASV ID to genome
    - asv_to_copies: Dict mapping ASV ID to number of copies (original sequences)
    """
    asv_to_genome = {}
    asv_to_copies = {}
    
    print(f"Reading deduplication mapping from {file_path}")
    df = pd.read_csv(file_path, sep='\t')
    
    for _, row in df.iterrows():
        asv_id = row['dedup_id']
        genome = row['genome']
        num_copies = row['num_original']
        
        asv_to_genome[asv_id] = genome
        asv_to_copies[asv_id] = num_copies
    
    print(f"  Loaded {len(asv_to_genome)} ASVs")
    print(f"  Example ASV IDs: {list(asv_to_genome.keys())[:3]}")
    
    return asv_to_genome, asv_to_copies


def read_genome_abundances(file_path: str) -> Dict[str, float]:
    """Read genome abundances from file."""
    print(f"Reading genome abundances from {file_path}")
    genome_abunds = {}
    
    df = pd.read_csv(file_path, sep='\t')
    for _, row in df.iterrows():
        genome = row['GenomeID']
        proportion = row['Proportion'] / 100.0  # Convert percentage to proportion
        genome_abunds[genome] = proportion
    
    print(f"  Loaded abundances for {len(genome_abunds)} genomes")
    
    return genome_abunds


def read_closest_file(file_path: str, max_distance: Optional[float] = None) -> Dict[str, int]:
    """
    Read closest matches file and count occurrences of each reference.
    Only counts reads with edit distance <= max_distance (if specified).
    """
    threshold_msg = f" with max distance {max_distance}" if max_distance is not None else ""
    print(f"Reading closest matches from {file_path}{threshold_msg}")
    
    ref_counts = Counter()
    total_reads = 0
    filtered_reads = 0
    
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                total_reads += 1
                try:
                    edit_distance = float(parts[1])
                    ref_id = parts[2]
                    
                    # Apply threshold if specified
                    if max_distance is None or edit_distance <= max_distance:
                        ref_counts[ref_id] += 1
                    else:
                        filtered_reads += 1
                        
                except ValueError:
                    continue  # Skip malformed lines
    
    if max_distance is not None:
        print(f"  Total reads: {total_reads}, Filtered out: {filtered_reads}, Kept: {total_reads - filtered_reads}")
    
    print(f"  Found counts for {len(ref_counts)} unique references")
    print(f"  Example reference IDs: {list(ref_counts.keys())[:3]}")
    
    return dict(ref_counts)


def create_simple_id_mapping(asv_to_genome: Dict[str, str]) -> Dict[str, str]:
    """
    Create mapping from original ASV IDs to simplified {genome}-{index} format.
    Since IDs are already in this format, just return identity mapping.
    """
    # The IDs are already simplified, so just return them as-is
    return {asv_id: asv_id for asv_id in asv_to_genome}


def calculate_expected_abundances(
    asv_to_genome: Dict[str, str],
    asv_to_copies: Dict[str, int],
    genome_abunds: Dict[str, float]
) -> Dict[str, float]:
    """
    Calculate expected abundances for each ASV.
    Formula: genome_proportion × (copy_count / total_copy_count_for_genome)
    """
    # Calculate total copies per genome
    genome_to_total_copies = defaultdict(int)
    for asv_id, genome in asv_to_genome.items():
        genome_to_total_copies[genome] += asv_to_copies[asv_id]
    
    # Calculate expected abundance for each ASV
    expected_abundances = {}
    for asv_id, genome in asv_to_genome.items():
        if genome in genome_abunds and genome_to_total_copies[genome] > 0:
            genome_abund = genome_abunds[genome]
            copy_count = asv_to_copies[asv_id]
            total_copies = genome_to_total_copies[genome]
            
            expected_abund = genome_abund * (copy_count / total_copies)
            expected_abundances[asv_id] = expected_abund
        else:
            # Genome not in abundance file, set to 0
            expected_abundances[asv_id] = 0
            if genome not in genome_abunds:
                print(f"  Warning: Genome '{genome}' not found in abundance file")
    
    print(f"  Calculated expected abundances for {len(expected_abundances)} ASVs")
    
    return expected_abundances


def get_normalized_counts(counts: Dict[str, int]) -> Dict[str, float]:
    """Convert raw counts to relative abundances (proportions)."""
    total_counts = sum(counts.values())
    if total_counts == 0:
        return {k: 0 for k in counts}
    
    return {k: v / total_counts for k, v in counts.items()}


def summarize_by_genome(abundances: Dict[str, float], asv_to_genome: Dict[str, str]) -> Dict[str, float]:
    """Summarize abundances by genome."""
    genome_abundances = defaultdict(float)
    
    for asv_id, abundance in abundances.items():
        if asv_id in asv_to_genome:
            genome = asv_to_genome[asv_id]
            genome_abundances[genome] += abundance
    
    return dict(genome_abundances)


def plot_expected_vs_observed(
    expected_abunds: Dict[str, float],
    real_abunds: Dict[str, float],
    sim_abunds: Dict[str, float],
    asv_to_genome: Dict[str, str],
    original_to_simple: Dict[str, str],
    out_dir: str,
    simplify_ids: bool = False,
    summarize_by_genome_flag: bool = False
):
    """Create plots comparing expected vs observed abundances."""
    os.makedirs(out_dir, exist_ok=True)
    
    # Set up publication-quality plot parameters
    plt.rcParams.update({
        'font.size': 14,
        'font.weight': 'bold',
        'axes.labelsize': 16,
        'axes.labelweight': 'bold',
        'axes.titlesize': 18,
        'axes.titleweight': 'bold',
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'legend.fontsize': 12,
        'figure.titlesize': 20,
        'figure.titleweight': 'bold'
    })
    
    # Summarize by genome if requested
    if summarize_by_genome_flag:
        expected_abunds = summarize_by_genome(expected_abunds, asv_to_genome)
        real_abunds = summarize_by_genome(real_abunds, asv_to_genome)
        sim_abunds = summarize_by_genome(sim_abunds, asv_to_genome)
        
        # For genome-level plot, we don't need the ASV-to-genome mapping
        plot_data = []
        for genome in sorted(expected_abunds.keys()):
            plot_data.append({
                'ID': genome,
                'Expected': expected_abunds.get(genome, 0),
                'Real Observed': real_abunds.get(genome, 0),
                'Simulated Observed': sim_abunds.get(genome, 0)
            })
    else:
        # Prepare ASV-level data for plotting
        plot_data = []
        for asv_id, expected in expected_abunds.items():
            genome = asv_to_genome[asv_id]
            
            # Display ID (original or simplified)
            display_id = original_to_simple[asv_id] if simplify_ids else asv_id
            
            # Get observed abundances (default to 0 if not found)
            real_observed = real_abunds.get(asv_id, 0)
            sim_observed = sim_abunds.get(asv_id, 0)
            
            # Add to plot data
            plot_data.append({
                'ID': display_id,
                'Genome': genome,
                'Expected': expected,
                'Real Observed': real_observed,
                'Simulated Observed': sim_observed
            })
    
    # Convert to DataFrame
    df = pd.DataFrame(plot_data)
    
    # Handle empty dataframe
    if len(df) == 0:
        print("Warning: No data to plot!")
        return
    
    # Sort appropriately
    if summarize_by_genome_flag:
        df = df.sort_values('Expected', ascending=False)
    else:
        df = df.sort_values(['Genome', 'Expected'], ascending=[True, False])
    
    # Calculate correlations (handle NaN)
    try:
        real_corr = np.corrcoef(df['Expected'], df['Real Observed'])[0, 1]
        if np.isnan(real_corr):
            real_corr = 0.0
    except:
        real_corr = 0.0
        
    try:
        sim_corr = np.corrcoef(df['Expected'], df['Simulated Observed'])[0, 1]
        if np.isnan(sim_corr):
            sim_corr = 0.0
    except:
        sim_corr = 0.0
    
    # 1. Combined bar plot
    fig_width = 20 if not summarize_by_genome_flag else 12
    plt.figure(figsize=(fig_width, 10))
    
    # Set positions for bars
    num_items = len(df)
    indices = np.arange(num_items)
    width = 0.25
    
    # Plot bars
    plt.bar(indices - width, df['Expected'], width, label='Expected', color='blue', alpha=0.8)
    plt.bar(indices, df['Real Observed'], width, label='Real Observed', color='green', alpha=0.8)
    plt.bar(indices + width, df['Simulated Observed'], width, label='Simulated Observed', color='red', alpha=0.8)
    
    # Add genome dividers and labels (only for ASV-level plots)
    if not summarize_by_genome_flag:
        genome_boundaries = []
        current_genome = None
        for i, genome in enumerate(df['Genome']):
            if genome != current_genome:
                if i > 0:
                    plt.axvline(x=i-0.5, color='black', linestyle='--', alpha=0.3)
                genome_boundaries.append((i, genome))
                current_genome = genome
        
        # Add genome labels at midpoints
        for i in range(len(genome_boundaries)):
            start_idx = genome_boundaries[i][0]
            genome = genome_boundaries[i][1]
            
            end_idx = num_items
            if i < len(genome_boundaries) - 1:
                end_idx = genome_boundaries[i+1][0]
            
            midpoint = (start_idx + end_idx - 1) / 2
            plt.text(midpoint, -0.02, genome, transform=plt.gca().get_xaxis_transform(),
                     ha='center', va='top', fontsize=12, fontweight='bold')
    
    # Label the x-axis
    plt.xticks(indices, df['ID'], rotation=90, fontsize=10)
    
    # Add correlations to the title
    level_text = "Genome-level" if summarize_by_genome_flag else "ASV-level"
    plt.title(f'{level_text} Expected vs Observed Abundances\nReal Correlation: {real_corr:.3f}, Sim Correlation: {sim_corr:.3f}')
    plt.ylabel('Relative Abundance')
    plt.xlabel('Genome' if summarize_by_genome_flag else 'ASV')
    plt.legend(frameon=True, fancybox=True, shadow=True)
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    
    suffix = '_genome' if summarize_by_genome_flag else ''
    plt.savefig(os.path.join(out_dir, f'expected_vs_observed_combined{suffix}.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Expected vs Real Observed Scatter Plot
    plt.figure(figsize=(10, 10))
    
    if summarize_by_genome_flag:
        # Simple scatter for genome-level
        plt.scatter(df['Expected'], df['Real Observed'], s=200, alpha=0.8, edgecolors='black', linewidth=1.5)
        
        # Add genome labels
        for _, row in df.iterrows():
            plt.annotate(row['ID'], (row['Expected'], row['Real Observed']), 
                        fontsize=10, ha='center', va='bottom')
    else:
        # Color by genome for ASV-level
        unique_genomes = df['Genome'].unique()
        genome_to_color = {genome: plt.cm.tab20(i % 20) for i, genome in enumerate(unique_genomes)}
        colors = [genome_to_color[genome] for genome in df['Genome']]
        
        plt.scatter(df['Expected'], df['Real Observed'], c=colors, alpha=0.8, s=100, edgecolors='black', linewidth=0.5)
        
        # Add legend for genomes
        handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10) 
                   for color in [genome_to_color[g] for g in unique_genomes]]
        plt.legend(handles, unique_genomes, title='Genome', loc='best', frameon=True)
    
    # Add reference line (y=x)
    max_val = max(df['Expected'].max(), df['Real Observed'].max()) * 1.1
    plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=2)
    
    plt.xlabel('Expected Abundance')
    plt.ylabel('Observed Abundance (Real)')
    plt.title(f'{level_text} Expected vs Real Observed\nCorrelation: {real_corr:.3f}')
    plt.grid(True, alpha=0.3)
    plt.xlim(-0.01, max_val)
    plt.ylim(-0.01, max_val)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'expected_vs_real_scatter{suffix}.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Expected vs Simulated Observed Scatter Plot
    plt.figure(figsize=(10, 10))
    
    if summarize_by_genome_flag:
        plt.scatter(df['Expected'], df['Simulated Observed'], s=200, alpha=0.8, edgecolors='black', linewidth=1.5)
        
        # Add genome labels
        for _, row in df.iterrows():
            plt.annotate(row['ID'], (row['Expected'], row['Simulated Observed']), 
                        fontsize=10, ha='center', va='bottom')
    else:
        plt.scatter(df['Expected'], df['Simulated Observed'], c=colors, alpha=0.8, s=100, edgecolors='black', linewidth=0.5)
        plt.legend(handles, unique_genomes, title='Genome', loc='best', frameon=True)
    
    # Add reference line (y=x)
    max_val = max(df['Expected'].max(), df['Simulated Observed'].max()) * 1.1
    plt.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=2)
    
    plt.xlabel('Expected Abundance')
    plt.ylabel('Observed Abundance (Simulated)')
    plt.title(f'{level_text} Expected vs Simulated Observed\nCorrelation: {sim_corr:.3f}')
    plt.grid(True, alpha=0.3)
    plt.xlim(-0.01, max_val)
    plt.ylim(-0.01, max_val)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f'expected_vs_sim_scatter{suffix}.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Save data to CSV
    df.to_csv(os.path.join(out_dir, f'expected_vs_observed_data{suffix}.csv'), index=False)
    
    print(f"Plots and data saved to {out_dir}")
    print(f"Correlation between Expected and Real Observed: {real_corr:.3f}")
    print(f"Correlation between Expected and Simulated Observed: {sim_corr:.3f}")


def main():
    args = parse_args()
    
    # Read input files
    asv_to_genome, asv_to_copies = read_dedup_mapping(args.dedup_mapping)
    genome_abunds = read_genome_abundances(args.genome_abunds)
    
    # Read observed counts using closest files with optional thresholds
    real_counts = read_closest_file(args.real_closest, args.real_closest_threshold)
    
    # For simulated data, use closest matches with optional threshold
    try:
        sim_counts = read_closest_file(args.sim_closest, args.sim_closest_threshold)
    except FileNotFoundError:
        print(f"Warning: Simulated data file not found: {args.sim_closest}")
        print("Using empty counts for simulated data")
        sim_counts = {}
    
    # Create mapping to simplified IDs if requested
    original_to_simple = create_simple_id_mapping(asv_to_genome)
    
    # Calculate expected abundances
    expected_abunds = calculate_expected_abundances(asv_to_genome, asv_to_copies, genome_abunds)
    
    # Normalize observed counts to relative abundances
    real_abunds = get_normalized_counts(real_counts)
    sim_abunds = get_normalized_counts(sim_counts)
    
    # Debug: Check if IDs match
    print("\nDebug: Checking ID matches...")
    print(f"Expected ASV IDs (first 5): {list(expected_abunds.keys())[:5]}")
    print(f"Real observed IDs (first 5): {list(real_abunds.keys())[:5]}")
    print(f"Sim observed IDs (first 5): {list(sim_abunds.keys())[:5]}")
    
    # Check for overlaps
    expected_ids = set(expected_abunds.keys())
    real_ids = set(real_abunds.keys())
    sim_ids = set(sim_abunds.keys())
    
    print(f"\nOverlap between expected and real: {len(expected_ids & real_ids)} / {len(expected_ids)}")
    print(f"Overlap between expected and sim: {len(expected_ids & sim_ids)} / {len(expected_ids)}")
    
    # Create plots
    plot_expected_vs_observed(
        expected_abunds, 
        real_abunds, 
        sim_abunds, 
        asv_to_genome,
        original_to_simple,
        args.out_dir,
        args.simplify_ids,
        args.summarize_by_genome
    )


if __name__ == "__main__":
    main()
