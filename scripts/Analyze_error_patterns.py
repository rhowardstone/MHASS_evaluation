#!/usr/bin/env python3
"""
Analyze_error_patterns.py
Analyze error patterns by examining how reads differ from their assigned ground truth.
For real data: reads vs consensus sequences
For simulated data: reads vs reference sequences
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import edlib
from collections import defaultdict, Counter
from tqdm import tqdm
import re
from multiprocessing import Pool, cpu_count
import itertools
import random

def parse_args():
    parser = argparse.ArgumentParser(description='Analyze error patterns in read assignments')
    parser.add_argument('--closest_file', required=True,
                        help='Closest file with read-to-reference assignments')
    parser.add_argument('--reads_fasta', required=True,
                        help='FASTA file with reads')
    parser.add_argument('--reference_fasta', required=True,
                        help='FASTA file with references (consensus for real, ref for sim)')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory for plots')
    parser.add_argument('--max_distance', type=int, default=None,
                        help='Maximum distance to analyze in detail (default: no limit)')
    parser.add_argument('--sample_size', type=int, default=None,
                        help='Number of reads to sample for analysis (default: all)')
    parser.add_argument('--threads', type=int, default=cpu_count(),
                        help='Number of threads to use (default: all CPUs)')
    parser.add_argument('--data_type', choices=['real', 'simulated'], required=True,
                        help='Type of data being analyzed')
    parser.add_argument('--random_seed', type=int, default=42,
                        help='Random seed for sampling (default: 42)')
    return parser.parse_args()

def read_fasta_dict(fasta_file):
    """Read FASTA file into dictionary, preserving full IDs."""
    sequences = {}
    with open(fasta_file, 'r') as f:
        current_id = None
        current_seq = []
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]  # Full ID including spaces
                current_seq = []
            elif line:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def analyze_single_alignment(args):
    """Analyze a single read-reference alignment. Used for parallel processing."""
    read_id, distance, ref_id, read_seq, ref_seq, max_distance = args
    
    # Skip if exceeds max distance (if specified)
    if max_distance is not None and distance > max_distance:
        return None, distance, None, None
    
    if not read_seq or not ref_seq:
        return None, distance, None, None
    
    # Analyze alignment
    result = edlib.align(read_seq, ref_seq, task="path", mode="NW")
    
    if result is None or 'cigar' not in result:
        return None, distance, None, None
    
    # Parse CIGAR string to count error types and positions
    cigar = result['cigar']
    subs, ins, dels = 0, 0, 0
    
    # Track errors by position
    errors_by_position = []
    ref_pos = 0
    
    # Extract all operations from CIGAR
    operations = re.findall(r'(\d+)([=XIDSHP])', cigar)
    
    for count_str, op_type in operations:
        count = int(count_str)
        if op_type == '=':  # Match
            ref_pos += count
        elif op_type == 'X':  # Mismatch (substitution)
            subs += count
            for i in range(count):
                errors_by_position.append((ref_pos + i, 'substitution'))
            ref_pos += count
        elif op_type == 'I':  # Insertion to target (deletion from query)
            dels += count
            for i in range(count):
                errors_by_position.append((ref_pos, 'deletion'))
        elif op_type == 'D':  # Deletion from target (insertion to query)
            ins += count
            for i in range(count):
                errors_by_position.append((ref_pos + i, 'insertion'))
            ref_pos += count
    
    return {'subs': subs, 'ins': ins, 'del': dels}, distance, distance, errors_by_position

def analyze_error_patterns_parallel(closest_file, reads_dict, refs_dict, max_distance=None, 
                                   sample_size=None, threads=1, random_seed=42):
    """Analyze error patterns from closest assignments using parallel processing."""
    
    # Read all entries from closest file
    all_work_items = []
    with open(closest_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            
            read_id = parts[0]
            distance = int(parts[1])
            ref_id = parts[2]
            
            # Get sequences
            read_seq = reads_dict.get(read_id)
            ref_seq = refs_dict.get(ref_id)
            
            if read_seq and ref_seq:
                all_work_items.append((read_id, distance, ref_id, read_seq, ref_seq, max_distance))
    
    # Random sampling if sample_size is specified
    if sample_size is not None and len(all_work_items) > sample_size:
        random.seed(random_seed)
        work_items = random.sample(all_work_items, sample_size)
        print(f"Randomly sampled {sample_size} reads from {len(all_work_items)} total reads")
    else:
        work_items = all_work_items
    
    # Process in parallel
    error_types = defaultdict(list)
    distance_distribution = []
    errors_by_distance = defaultdict(lambda: {'subs': 0, 'ins': 0, 'del': 0, 'total': 0})
    position_errors = defaultdict(lambda: {'subs': 0, 'ins': 0, 'del': 0})
    
    with Pool(threads) as pool:
        # Process with progress bar
        results = list(tqdm(
            pool.imap(analyze_single_alignment, work_items, chunksize=100),
            total=len(work_items),
            desc="Analyzing alignments"
        ))
    
    # Aggregate results
    for error_data, distance, dist_for_analysis, errors_by_pos in results:
        if distance is not None:
            distance_distribution.append(distance)
        
        if error_data is not None and dist_for_analysis is not None:
            error_types['substitutions'].append(error_data['subs'])
            error_types['insertions'].append(error_data['ins'])
            error_types['deletions'].append(error_data['del'])
            
            errors_by_distance[dist_for_analysis]['subs'] += error_data['subs']
            errors_by_distance[dist_for_analysis]['ins'] += error_data['ins']
            errors_by_distance[dist_for_analysis]['del'] += error_data['del']
            errors_by_distance[dist_for_analysis]['total'] += 1
            
            # Aggregate position-based errors
            if errors_by_pos:
                for pos, error_type in errors_by_pos:
                    if error_type == 'substitution':
                        position_errors[pos]['subs'] += 1
                    elif error_type == 'insertion':
                        position_errors[pos]['ins'] += 1
                    elif error_type == 'deletion':
                        position_errors[pos]['del'] += 1
    
    return error_types, distance_distribution, errors_by_distance, position_errors

def plot_position_based_errors(position_errors, output_dir, data_type):
    """Create plots showing error patterns by position in alignment."""
    if not position_errors:
        print("No position-based error data to plot")
        return
    
    # Convert to sorted lists
    positions = sorted(position_errors.keys())
    subs_counts = [position_errors[pos]['subs'] for pos in positions]
    ins_counts = [position_errors[pos]['ins'] for pos in positions]
    del_counts = [position_errors[pos]['del'] for pos in positions]
    
    # 1. Raw counts stacked plot
    plt.figure(figsize=(14, 8))
    plt.stackplot(positions, subs_counts, ins_counts, del_counts,
                  labels=['Substitutions', 'Insertions', 'Deletions'],
                  colors=['#D62728', '#2CA02C', '#1F77B4'],  # Red, Green, Blue - high contrast
                  alpha=0.8)
    plt.xlabel('Position in Alignment (bp)', fontsize=12)
    plt.ylabel('Number of Errors', fontsize=12)
    plt.title(f'Error Distribution by Position (Raw Counts) - {data_type} data', fontsize=14)
    plt.legend(loc='upper right')
    plt.grid(True, alpha=0.3)
    
    # Add annotation if there's a notable peak around 200bp
    if 180 <= max(positions, key=lambda x: sum([position_errors[x][k] for k in ['subs', 'ins', 'del']])) <= 220:
        peak_pos = max(positions, key=lambda x: sum([position_errors[x][k] for k in ['subs', 'ins', 'del']]))
        peak_height = sum([position_errors[peak_pos][k] for k in ['subs', 'ins', 'del']])
        plt.annotate(f'Peak at {peak_pos}bp', xy=(peak_pos, peak_height), 
                    xytext=(peak_pos + 50, peak_height + 10),
                    arrowprops=dict(arrowstyle='->', color='red'),
                    fontsize=10, color='red')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'errors_by_position_raw_{data_type}.png'), dpi=300)
    plt.close()
    
    # 2. Proportions stacked plot
    plt.figure(figsize=(14, 8))
    
    # Calculate proportions
    total_errors = [subs_counts[i] + ins_counts[i] + del_counts[i] for i in range(len(positions))]
    subs_prop = []
    ins_prop = []
    del_prop = []
    
    for i in range(len(positions)):
        if total_errors[i] > 0:
            subs_prop.append(subs_counts[i] / total_errors[i])
            ins_prop.append(ins_counts[i] / total_errors[i])
            del_prop.append(del_counts[i] / total_errors[i])
        else:
            subs_prop.append(0)
            ins_prop.append(0)
            del_prop.append(0)
    
    plt.stackplot(positions, subs_prop, ins_prop, del_prop,
                  labels=['Substitutions', 'Insertions', 'Deletions'],
                  colors=['#D62728', '#2CA02C', '#1F77B4'],  # Red, Green, Blue - high contrast
                  alpha=0.8)
    plt.xlabel('Position in Alignment (bp)', fontsize=12)
    plt.ylabel('Proportion of Error Type', fontsize=12)
    plt.title(f'Error Composition by Position (Proportions) - {data_type} data', fontsize=14)
    plt.legend(loc='upper right')
    plt.ylim(0, 1)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'errors_by_position_prop_{data_type}.png'), dpi=300)
    plt.close()

def plot_error_analysis(error_types, distance_dist, errors_by_dist, output_dir, data_type):
    """Create plots showing error patterns."""
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Error type distribution - using violin plot instead of boxplot
    plt.figure(figsize=(10, 6))
    data = []
    for error_type, counts in error_types.items():
        for count in counts:
            data.append({'Error Type': error_type.capitalize(), 'Count': count})
    
    df = pd.DataFrame(data)
    if len(df) > 0:
        # Use consistent colors for error types
        palette = {'Substitutions': '#D62728', 'Insertions': '#2CA02C', 'Deletions': '#1F77B4'}
        sns.violinplot(data=df, x='Error Type', y='Count', inner='quartile', palette=palette)
        plt.title(f'Distribution of Error Types - {data_type} data')
        plt.ylabel('Number of Errors per Read')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'error_types_{data_type}.png'), dpi=300)
    plt.close()
    
    # 2. Error composition by distance - no x-axis limit
    distances = sorted(errors_by_dist.keys())
    subs_prop = []
    ins_prop = []
    del_prop = []
    
    for d in distances:
        total_errors = (errors_by_dist[d]['subs'] + 
                       errors_by_dist[d]['ins'] + 
                       errors_by_dist[d]['del'])
        if total_errors > 0:
            subs_prop.append(errors_by_dist[d]['subs'] / total_errors)
            ins_prop.append(errors_by_dist[d]['ins'] / total_errors)
            del_prop.append(errors_by_dist[d]['del'] / total_errors)
        else:
            subs_prop.append(0)
            ins_prop.append(0)
            del_prop.append(0)
    
    if distances:
        plt.figure(figsize=(12, 6))
        plt.stackplot(distances, subs_prop, ins_prop, del_prop, 
                      labels=['Substitutions', 'Insertions', 'Deletions'],
                      colors=['#D62728', '#2CA02C', '#1F77B4'],  # Red, Green, Blue - high contrast
                      alpha=0.8)
        plt.xlabel('Edit Distance')
        plt.ylabel('Proportion of Error Type')
        plt.title(f'Error Composition by Distance - {data_type} data')
        plt.legend(loc='upper left')
        # No plt.xlim() - show full range
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'error_composition_by_distance_{data_type}.png'), dpi=300)
        plt.close()
    
    # 3. Distance distribution - no x-axis limit
    if distance_dist:
        plt.figure(figsize=(10, 6))
        plt.hist(distance_dist, bins=50, edgecolor='black', alpha=0.7)
        plt.xlabel('Edit Distance')
        plt.ylabel('Number of Reads')
        plt.title(f'Distribution of Edit Distances - {data_type} data')
        # No plt.xlim() - show full range
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'distance_distribution_{data_type}.png'), dpi=300)
        plt.close()

def main():
    args = parse_args()
    
    print(f"Loading sequences from {args.reads_fasta}")
    reads_dict = read_fasta_dict(args.reads_fasta)
    print(f"Loaded {len(reads_dict)} reads")
    
    print(f"Loading references from {args.reference_fasta}")
    refs_dict = read_fasta_dict(args.reference_fasta)
    print(f"Loaded {len(refs_dict)} references")
    
    print(f"Analyzing error patterns using {args.threads} threads...")
    if args.sample_size:
        print(f"Using random sampling with {args.sample_size} reads (seed: {args.random_seed})")
    if args.max_distance:
        print(f"Only analyzing reads with distance <= {args.max_distance}")
    
    error_types, distance_dist, errors_by_dist, position_errors = analyze_error_patterns_parallel(
        args.closest_file, reads_dict, refs_dict, 
        args.max_distance, args.sample_size, args.threads, args.random_seed
    )
    
    print("Creating plots...")
    plot_error_analysis(error_types, distance_dist, errors_by_dist, 
                       args.output_dir, args.data_type)
    
    # Create position-based error plots
    print("Creating position-based error plots...")
    plot_position_based_errors(position_errors, args.output_dir, args.data_type)
    
    # Print summary
    print(f"\nError Summary for {args.data_type} data:")
    print(f"Total reads analyzed: {len(distance_dist)}")
    if distance_dist:
        print(f"Mean distance: {np.mean(distance_dist):.1f}")
        print(f"Median distance: {np.median(distance_dist):.1f}")
        print(f"Max distance: {max(distance_dist)}")
    
    total_subs = sum(error_types['substitutions'])
    total_ins = sum(error_types['insertions'])
    total_dels = sum(error_types['deletions'])
    total_errors = total_subs + total_ins + total_dels
    
    if total_errors > 0:
        print(f"\nError composition:")
        print(f"  Substitutions: {total_subs} ({total_subs/total_errors*100:.1f}%)")
        print(f"  Insertions: {total_ins} ({total_ins/total_errors*100:.1f}%)")
        print(f"  Deletions: {total_dels} ({total_dels/total_errors*100:.1f}%)")
    
    # Position-based summary
    if position_errors:
        max_error_pos = max(position_errors.keys(), 
                           key=lambda x: sum(position_errors[x].values()))
        max_errors = sum(position_errors[max_error_pos].values())
        print(f"\nPosition with most errors: {max_error_pos}bp ({max_errors} errors)")
        
        # Check for peaks around 200bp
        peak_region = {pos: sum(position_errors[pos].values()) 
                      for pos in position_errors if 180 <= pos <= 220}
        if peak_region:
            peak_pos = max(peak_region, key=peak_region.get)
            print(f"Peak in 180-220bp region: {peak_pos}bp ({peak_region[peak_pos]} errors)")

if __name__ == '__main__':
    main()
