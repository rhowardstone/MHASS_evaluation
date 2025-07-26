#!/usr/bin/env python3
"""
Check_ASV_misassignments.py
Analyze misassignments at the ASV level and create confusion matrix showing
which ASVs are getting confused with each other.
"""

import argparse
import os
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Analyze ASV-level misassignments')
    parser.add_argument('--closest_file', required=True,
                        help='Closest file (e.g., ref-closest-sim.txt or consensus-closest-real.txt)')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory for plots and reports')
    parser.add_argument('--max_distance', type=int, default=50,
                        help='Maximum edit distance to consider (default: 50)')
    parser.add_argument('--data_type', choices=['real', 'simulated'], required=True,
                        help='Type of data being analyzed')
    parser.add_argument('--min_reads', type=int, default=5,
                        help='Minimum reads per ASV to include in confusion matrix (default: 5)')
    parser.add_argument('--sequence_file', help='Optional: simulated reads file to extract true source from read IDs')
    parser.add_argument('--top_asvs', type=int, default=50,
                        help='Number of top ASVs to include in confusion matrix (default: 50)')
    return parser.parse_args()

def extract_genome(seq_id):
    """Extract genome from deduped ID format like 'Escherichia_coli-1-2-3'."""
    if pd.isna(seq_id) or not isinstance(seq_id, str):
        return None
    
    # Find the position of the first dash followed by a digit
    match = re.search(r'^(.+?)-\d', seq_id)
    if match:
        return match.group(1)
    return seq_id

def extract_template_asv_from_read_id(read_id, template_to_asv_map=None):
    """
    Extract the true source ASV from a simulated read ID.
    Expected format: something like "template_XXX_YYY" where template_XXX maps to an ASV
    """
    # Try to extract template ID from read name
    template_match = re.search(r'template_(\d+)', read_id)
    if template_match and template_to_asv_map:
        template_id = f"template_{template_match.group(1)}"
        return template_to_asv_map.get(template_id, None)
    
    # Alternative: if read ID contains the ASV ID directly
    # Check if read ID contains any known ASV pattern
    asv_match = re.search(r'([A-Z][a-z]+_[a-z]+(?:-\d+)+)', read_id)
    if asv_match:
        return asv_match.group(1)
    
    return None

def build_template_to_asv_map(sequence_file_mapping):
    """
    Build a mapping from template IDs to ASV IDs using the sequence file mapping.
    This would need the actual mapping file from the simulation.
    """
    # For now, return empty map - this would need to be implemented based on
    # how your simulation tracks which templates come from which ASVs
    return {}

def load_closest_with_truth(filename, max_distance=None, template_to_asv_map=None):
    """Load closest file and extract true source ASV for simulated data."""
    data = []
    
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
                
            read_id = parts[0]
            edit_distance = int(parts[1])
            
            # Skip if distance exceeds threshold
            if max_distance is not None and edit_distance > max_distance:
                continue
            
            # Try to extract true source ASV from read ID (for simulated data)
            true_asv = extract_template_asv_from_read_id(read_id, template_to_asv_map)
            
            # Get assigned reference(s) - we'll take the first one
            assigned_asv = parts[2]
            
            data.append({
                'read_id': read_id,
                'true_asv': true_asv,
                'assigned_asv': assigned_asv,
                'assigned_genome': extract_genome(assigned_asv),
                'edit_distance': edit_distance
            })
    
    return pd.DataFrame(data)

def analyze_asv_confusion(df):
    """Analyze confusion between ASVs."""
    # For simulated data with known truth
    if 'true_asv' in df.columns and df['true_asv'].notna().any():
        # Build confusion matrix: true ASV vs assigned ASV
        confusion = defaultdict(lambda: defaultdict(int))
        
        for _, row in df.iterrows():
            if pd.notna(row['true_asv']):
                confusion[row['true_asv']][row['assigned_asv']] += 1
        
        return confusion, 'truth_based'
    
    # For real data or when truth is unknown, use distance-based analysis
    # Look at pairs of ASVs that have very similar sequences (low edit distances)
    confusion = defaultdict(lambda: defaultdict(int))
    
    # Group by assigned ASV and look at distance distributions
    for asv in df['assigned_asv'].unique():
        asv_data = df[df['assigned_asv'] == asv]
        
        # Count assignments by edit distance
        for _, row in asv_data.iterrows():
            # For very low edit distances, increment confusion count
            if row['edit_distance'] < 10:  # Threshold for "could be confused"
                confusion[asv][asv] += 1
    
    return confusion, 'distance_based'

def plot_asv_confusion_matrix(confusion_dict, output_dir, data_type, top_n=50):
    """Create a heatmap showing ASV-to-ASV confusion matrix."""
    # Convert to DataFrame
    all_asvs = sorted(set(list(confusion_dict.keys()) + 
                          [asv for targets in confusion_dict.values() for asv in targets.keys()]))
    
    # Select top ASVs by total reads
    asv_totals = {}
    for source_asv in all_asvs:
        total = sum(confusion_dict.get(source_asv, {}).values())
        # Also add reads assigned to this ASV from other sources
        for targets in confusion_dict.values():
            total += targets.get(source_asv, 0)
        asv_totals[source_asv] = total
    
    # Get top N ASVs
    top_asvs = sorted(asv_totals.keys(), key=lambda x: asv_totals[x], reverse=True)[:top_n]
    
    # Build matrix for top ASVs
    matrix_data = []
    for source in top_asvs:
        row = []
        for target in top_asvs:
            count = confusion_dict.get(source, {}).get(target, 0)
            row.append(count)
        matrix_data.append(row)
    
    confusion_df = pd.DataFrame(matrix_data, index=top_asvs, columns=top_asvs)
    
    # Calculate accuracy (diagonal vs off-diagonal)
    diagonal_sum = np.trace(confusion_df.values)
    total_sum = confusion_df.values.sum()
    accuracy = diagonal_sum / total_sum if total_sum > 0 else 0
    
    # Create heatmap
    plt.figure(figsize=(16, 14))
    
    # Use log scale for better visualization
    # Add 1 to avoid log(0)
    log_confusion = np.log10(confusion_df + 1)
    
    sns.heatmap(
        log_confusion,
        annot=confusion_df.values,  # Show actual counts
        fmt='d',
        cmap='YlOrRd',
        cbar_kws={'label': 'log10(Reads + 1)'},
        square=True
    )
    
    plt.title(f'ASV Confusion Matrix - {data_type} data\n' + 
              f'Top {top_n} ASVs by read count (Accuracy: {accuracy:.2%})')
    plt.xlabel('Assigned ASV')
    plt.ylabel('True ASV' if confusion_dict else 'Source ASV')
    
    # Rotate labels for better readability
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, f'asv_confusion_matrix_{data_type}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Also save the confusion matrix as CSV
    csv_file = os.path.join(output_dir, f'asv_confusion_matrix_{data_type}.csv')
    confusion_df.to_csv(csv_file)
    
    return confusion_df, accuracy

def plot_misassignment_patterns(df, confusion_df, output_dir, data_type):
    """Create additional plots showing misassignment patterns."""
    # 1. Bar plot of most confused ASV pairs
    off_diagonal = []
    
    for i, source in enumerate(confusion_df.index):
        for j, target in enumerate(confusion_df.columns):
            if i != j and confusion_df.iloc[i, j] > 0:
                off_diagonal.append({
                    'pair': f'{source} → {target}',
                    'count': confusion_df.iloc[i, j],
                    'source_genome': extract_genome(source),
                    'target_genome': extract_genome(target),
                    'same_genome': extract_genome(source) == extract_genome(target)
                })
    
    if off_diagonal:
        off_diag_df = pd.DataFrame(off_diagonal)
        off_diag_df = off_diag_df.sort_values('count', ascending=False).head(30)
        
        plt.figure(figsize=(12, 8))
        colors = ['blue' if x else 'red' for x in off_diag_df['same_genome']]
        
        plt.barh(range(len(off_diag_df)), off_diag_df['count'], color=colors)
        plt.yticks(range(len(off_diag_df)), off_diag_df['pair'])
        plt.xlabel('Number of Misassigned Reads')
        plt.title(f'Top Misassigned ASV Pairs - {data_type} data\n' +
                  'Blue: Same genome, Red: Different genomes')
        plt.tight_layout()
        
        output_file = os.path.join(output_dir, f'top_misassignments_{data_type}.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
    
    # 2. Summary plot by genome
    genome_confusion = defaultdict(lambda: defaultdict(int))
    
    for source in confusion_df.index:
        source_genome = extract_genome(source)
        for target in confusion_df.columns:
            target_genome = extract_genome(target)
            count = confusion_df.loc[source, target]
            genome_confusion[source_genome][target_genome] += count
    
    # Convert to DataFrame
    genomes = sorted(set(list(genome_confusion.keys()) + 
                        [g for targets in genome_confusion.values() for g in targets.keys()]))
    
    genome_matrix = []
    for source in genomes:
        row = []
        for target in genomes:
            count = genome_confusion[source][target]
            row.append(count)
        genome_matrix.append(row)
    
    genome_df = pd.DataFrame(genome_matrix, index=genomes, columns=genomes)
    
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        genome_df,
        annot=True,
        fmt='d',
        cmap='Blues',
        cbar_kws={'label': 'Total Reads'},
        square=True
    )
    plt.title(f'Genome-level Assignment Summary - {data_type} data')
    plt.xlabel('Assigned Genome')
    plt.ylabel('True/Source Genome')
    plt.tight_layout()
    
    output_file = os.path.join(output_dir, f'genome_summary_{data_type}.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def generate_detailed_report(df, confusion_df, accuracy, output_dir, data_type):
    """Generate a detailed report of misassignments."""
    report_file = os.path.join(output_dir, f'misassignment_report_{data_type}.txt')
    
    with open(report_file, 'w') as f:
        f.write(f"ASV Misassignment Analysis Report - {data_type} data\n")
        f.write("=" * 60 + "\n\n")
        
        f.write(f"Total assignments analyzed: {len(df)}\n")
        f.write(f"Unique ASVs: {df['assigned_asv'].nunique()}\n")
        f.write(f"Overall accuracy: {accuracy:.2%}\n\n")
        
        # Find most problematic ASVs
        f.write("Most Misassigned ASVs:\n")
        f.write("-" * 40 + "\n")
        
        misassignment_rates = {}
        for asv in confusion_df.index:
            total = confusion_df.loc[asv].sum()
            correct = confusion_df.loc[asv, asv]
            if total > 0:
                error_rate = 1 - (correct / total)
                misassignment_rates[asv] = {
                    'error_rate': error_rate,
                    'total': total,
                    'correct': correct,
                    'incorrect': total - correct
                }
        
        # Sort by error rate
        sorted_asvs = sorted(misassignment_rates.items(), 
                            key=lambda x: x[1]['error_rate'], 
                            reverse=True)
        
        for asv, stats in sorted_asvs[:20]:
            if stats['error_rate'] > 0:
                f.write(f"\n{asv} ({extract_genome(asv)}):\n")
                f.write(f"  Error rate: {stats['error_rate']:.1%}\n")
                f.write(f"  Total reads: {stats['total']}\n")
                f.write(f"  Misassigned: {stats['incorrect']}\n")
                
                # Show where they went
                f.write("  Misassigned to:\n")
                for target in confusion_df.columns:
                    if target != asv:
                        count = confusion_df.loc[asv, target]
                        if count > 0:
                            target_genome = extract_genome(target)
                            f.write(f"    {target} ({target_genome}): {count} reads\n")
        
        # Cross-genome misassignments
        f.write("\n\nCross-Genome Misassignments:\n")
        f.write("-" * 40 + "\n")
        
        cross_genome = []
        for source in confusion_df.index:
            source_genome = extract_genome(source)
            for target in confusion_df.columns:
                if source != target:
                    target_genome = extract_genome(target)
                    count = confusion_df.loc[source, target]
                    if count > 0 and source_genome != target_genome:
                        cross_genome.append({
                            'source': source,
                            'target': target,
                            'count': count,
                            'source_genome': source_genome,
                            'target_genome': target_genome
                        })
        
        cross_genome.sort(key=lambda x: x['count'], reverse=True)
        
        for item in cross_genome[:20]:
            f.write(f"\n{item['source']} ({item['source_genome']}) → " +
                   f"{item['target']} ({item['target_genome']}): {item['count']} reads\n")

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"Loading closest file: {args.closest_file}")
    
    # For simulated data, try to build template-to-ASV mapping
    template_to_asv_map = {}
    if args.data_type == 'simulated' and args.sequence_file:
        template_to_asv_map = build_template_to_asv_map(args.sequence_file)
    
    df = load_closest_with_truth(args.closest_file, args.max_distance, template_to_asv_map)
    print(f"Loaded {len(df)} assignments with distance <= {args.max_distance}")
    
    # Analyze ASV confusion
    print("Analyzing ASV-level confusion...")
    confusion_dict, analysis_type = analyze_asv_confusion(df)
    
    # Create confusion matrix plot
    print("Creating ASV confusion matrix...")
    confusion_df, accuracy = plot_asv_confusion_matrix(
        confusion_dict, args.output_dir, args.data_type, args.top_asvs
    )
    
    # Create additional analysis plots
    print("Creating misassignment pattern plots...")
    plot_misassignment_patterns(df, confusion_df, args.output_dir, args.data_type)
    
    # Generate detailed report
    print("Generating detailed report...")
    generate_detailed_report(df, confusion_df, accuracy, args.output_dir, args.data_type)
    
    print(f"\nAnalysis complete! Results saved to {args.output_dir}")
    print(f"Overall accuracy: {accuracy:.2%}")
    
    # Print summary of worst offenders
    print("\nTop 5 most confused ASV pairs:")
    off_diagonal = []
    for i, source in enumerate(confusion_df.index):
        for j, target in enumerate(confusion_df.columns):
            if i != j and confusion_df.iloc[i, j] > 0:
                off_diagonal.append((source, target, confusion_df.iloc[i, j]))
    
    off_diagonal.sort(key=lambda x: x[2], reverse=True)
    for source, target, count in off_diagonal[:5]:
        source_genome = extract_genome(source)
        target_genome = extract_genome(target)
        same = "same genome" if source_genome == target_genome else "DIFFERENT GENOMES"
        print(f"  {source} → {target}: {count} reads ({same})")

if __name__ == '__main__':
    main()
