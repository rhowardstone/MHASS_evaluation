#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.auto import tqdm
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description='Compare sequence length distributions between two FASTA files.')
    parser.add_argument('-1', '--fasta1', required=True,
                        help='First FASTA file (e.g., reference sequences)')
    parser.add_argument('-2', '--fasta2', required=True,
                        help='Second FASTA file (e.g., consensus sequences)')
    parser.add_argument('-n1', '--name1', default='Reference',
                        help='Label for the first dataset (default: Reference)')
    parser.add_argument('-n2', '--name2', default='Consensus',
                        help='Label for the second dataset (default: Consensus)')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Output directory for plots')
    parser.add_argument('-p', '--prefix', default='length_distribution',
                        help='Prefix for output files (default: length_distribution)')
    parser.add_argument('--bins', type=int, default=50,
                        help='Number of bins for histograms (default: 50)')
    parser.add_argument('--min_length', type=int, default=None,
                        help='Minimum sequence length to include (default: None)')
    parser.add_argument('--max_length', type=int, default=None,
                        help='Maximum sequence length to include (default: None)')
    return parser.parse_args()

def read_fasta_lengths(fasta_file, min_length=None, max_length=None):
    """Read a FASTA file and return a dictionary of sequence IDs and their lengths."""
    print(f"Reading sequences from {fasta_file}")
    lengths = {}
    
    for record in tqdm(SeqIO.parse(fasta_file, "fasta"), desc="Reading sequences"):
        seq_length = len(record.seq)
        
        # Apply length filters if specified
        if min_length is not None and seq_length < min_length:
            continue
        if max_length is not None and seq_length > max_length:
            continue
            
        lengths[record.id] = seq_length
    
    print(f"Found {len(lengths)} sequences in {fasta_file}")
    return lengths

def analyze_length_distribution(lengths):
    """Calculate statistics for sequence lengths."""
    if not lengths:
        return None
        
    values = list(lengths.values())
    stats = {
        'count': len(values),
        'mean': np.mean(values),
        'median': np.median(values),
        'std': np.std(values),
        'min': min(values),
        'max': max(values),
        'range': max(values) - min(values),
        'q1': np.percentile(values, 25),
        'q3': np.percentile(values, 75),
        'iqr': np.percentile(values, 75) - np.percentile(values, 25)
    }
    return stats

def plot_length_distributions(lengths1, lengths2, name1, name2, output_dir, prefix, bins=50):
    """Create plots comparing sequence length distributions."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Convert to DataFrames for easier plotting
    df1 = pd.DataFrame({'length': list(lengths1.values()), 'source': name1})
    df2 = pd.DataFrame({'length': list(lengths2.values()), 'source': name2})
    combined = pd.concat([df1, df2])
    
    # Calculate statistics
    stats1 = analyze_length_distribution(lengths1)
    stats2 = analyze_length_distribution(lengths2)
    
    if not stats1 or not stats2:
        print("Error: One or both datasets are empty after filtering.")
        return
    
    print(f"\n{name1} Length Statistics:")
    for key, value in stats1.items():
        print(f"  {key}: {value:.2f}")
    
    print(f"\n{name2} Length Statistics:")
    for key, value in stats2.items():
        print(f"  {key}: {value:.2f}")
    
    # Save statistics to CSV
    stats_df = pd.DataFrame({
        name1: [stats1['count'], stats1['mean'], stats1['median'], stats1['std'], 
                stats1['min'], stats1['max'], stats1['range'], 
                stats1['q1'], stats1['q3'], stats1['iqr']],
        name2: [stats2['count'], stats2['mean'], stats2['median'], stats2['std'], 
                stats2['min'], stats2['max'], stats2['range'], 
                stats2['q1'], stats2['q3'], stats2['iqr']]
    }, index=['Count', 'Mean', 'Median', 'Std Dev', 'Min', 'Max', 'Range', 
              'Q1 (25%)', 'Q3 (75%)', 'IQR'])
    
    stats_df.to_csv(os.path.join(output_dir, f"{prefix}_statistics.csv"))
    
    # 1. Histogram with density curve (separate)
    plt.figure(figsize=(10, 6))
    sns.histplot(df1['length'], kde=True, color='blue', label=name1)
    plt.xlabel("Sequence Length (bp)")
    plt.ylabel("Count")
    plt.title(f"Sequence Length Distribution: {name1}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{prefix}_{name1.lower()}_histogram.png"), dpi=300)
    plt.close()
    
    plt.figure(figsize=(10, 6))
    sns.histplot(df2['length'], kde=True, color='red', label=name2)
    plt.xlabel("Sequence Length (bp)")
    plt.ylabel("Count")
    plt.title(f"Sequence Length Distribution: {name2}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{prefix}_{name2.lower()}_histogram.png"), dpi=300)
    plt.close()
    
    # 2. Overlaid density plots
    plt.figure(figsize=(12, 6))
    sns.kdeplot(df1['length'], color='blue', label=name1, fill=True, alpha=0.3)
    sns.kdeplot(df2['length'], color='red', label=name2, fill=True, alpha=0.3)
    plt.xlabel("Sequence Length (bp)")
    plt.ylabel("Density")
    plt.title(f"Sequence Length Distribution: {name1} vs {name2}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{prefix}_density_overlay.png"), dpi=300)
    plt.close()
    
    # 3. Overlaid histograms with transparency
    plt.figure(figsize=(12, 6))
    sns.histplot(df1['length'], bins=bins, color='blue', alpha=0.5, label=name1, kde=True)
    sns.histplot(df2['length'], bins=bins, color='red', alpha=0.5, label=name2, kde=True)
    plt.xlabel("Sequence Length (bp)")
    plt.ylabel("Count")
    plt.title(f"Sequence Length Distribution: {name1} vs {name2}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{prefix}_histogram_overlay.png"), dpi=300)
    plt.close()
    
    # 4. Box plot comparison
    plt.figure(figsize=(8, 6))
    sns.boxplot(x='source', y='length', data=combined, palette=['blue', 'red'])
    plt.xlabel("Dataset")
    plt.ylabel("Sequence Length (bp)")
    plt.title(f"Sequence Length Comparison: {name1} vs {name2}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{prefix}_boxplot.png"), dpi=300)
    plt.close()
    
    # 5. Violin plot comparison
    plt.figure(figsize=(8, 6))
    sns.violinplot(x='source', y='length', data=combined, palette=['blue', 'red'], inner='quartile')
    plt.xlabel("Dataset")
    plt.ylabel("Sequence Length (bp)")
    plt.title(f"Sequence Length Comparison: {name1} vs {name2}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{prefix}_violinplot.png"), dpi=300)
    plt.close()
    
    # 6. CDF plot
    plt.figure(figsize=(10, 6))
    for data, label, color in zip([df1['length'], df2['length']], [name1, name2], ['blue', 'red']):
        sorted_data = np.sort(data)
        y = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
        plt.plot(sorted_data, y, label=label, color=color)
    
    plt.xlabel("Sequence Length (bp)")
    plt.ylabel("Cumulative Proportion")
    plt.title(f"Cumulative Distribution: {name1} vs {name2}")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{prefix}_cdf.png"), dpi=300)
    plt.close()
    
    print(f"Plots saved to {output_dir}")

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Read sequence lengths from FASTA files
    lengths1 = read_fasta_lengths(args.fasta1, args.min_length, args.max_length)
    lengths2 = read_fasta_lengths(args.fasta2, args.min_length, args.max_length)
    
    # Plot length distributions
    plot_length_distributions(
        lengths1, 
        lengths2, 
        args.name1, 
        args.name2, 
        args.output_dir, 
        args.prefix,
        args.bins
    )
    
    print("Analysis complete!")

if __name__ == '__main__':
    main()
