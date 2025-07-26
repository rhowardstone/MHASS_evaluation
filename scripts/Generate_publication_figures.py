#!/usr/bin/env python3
"""
generate_publication_figures.py
Generate all publication figures for MHASS paper with consistent formatting
Updated to generate error position data on the fly with full parallel processing
"""

import matplotlib
matplotlib.use('Agg')  # Set backend before importing pyplot

import argparse
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from matplotlib.gridspec import GridSpec
import re
import json
from scipy import stats
from collections import Counter, defaultdict
from pathlib import Path
import edlib
from Bio import SeqIO
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import random

# Set global font properties for bold text
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11
plt.rcParams['font.size'] = 12

def parse_args():
    parser = argparse.ArgumentParser(description='Generate all publication figures for MHASS')
    parser.add_argument('--output_dir', default='./Publication_figures',
                        help='Output directory for figures')
    parser.add_argument('--dpi', type=int, default=600,
                        help='DPI for output figures (default: 600)')
    parser.add_argument('--threads', type=int, default=None,
                        help='Number of threads for parallel processing (default: all CPUs)')
    return parser.parse_args()

def extract_genome(seq_id):
    """Extract genome from deduped ID format like 'Escherichia_coli-1-2-3'."""
    if pd.isna(seq_id) or not isinstance(seq_id, str):
        return None
    
    match = re.search(r'^(.+?)-\d', seq_id)
    if match:
        return match.group(1)
    return seq_id

def format_genome_name(genome_name):
    """Format genome name to species abbreviation format (e.g., E. coli) with italics"""
    if pd.isna(genome_name) or not isinstance(genome_name, str):
        return genome_name
    
    name = genome_name.replace('_', ' ')
    parts = name.split()
    
    if len(parts) >= 2:
        genus = parts[0]
        species = ' '.join(parts[1:])
        # Return with italic formatting for matplotlib
        return f"$\\it{{{genus[0]}.\\ {species}}}$"
    else:
        return f"$\\it{{{name}}}$"

def read_closest(closest_file, max_distance=None):
    """Read the closest.txt file to get edit distances."""
    data = []
    with open(closest_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            read_id = parts[0]
            edit_distance = int(parts[1])
            
            if max_distance is not None and edit_distance > max_distance:
                continue
                
            closest_refs = parts[2:]
            for ref in closest_refs:
                data.append((read_id, ref, edit_distance))
    return pd.DataFrame(data, columns=['read_id', 'reference_id', 'edit_distance'])

def read_fasta_dict(fasta_file):
    """Read FASTA file into dictionary."""
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
    """Analyze a single read-reference alignment for error patterns."""
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

def generate_error_position_data_parallel(closest_file, reads_fasta, reference_fasta, data_type, threads=None):
    """Generate error position data on the fly with parallel processing."""
    print(f"    Generating error position data for {data_type}...")
    
    # Set thread count
    if threads is None:
        threads = cpu_count()
    
    # Read sequences
    print(f"      Reading sequences...")
    reads_dict = read_fasta_dict(reads_fasta)
    refs_dict = read_fasta_dict(reference_fasta)
    
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
                all_work_items.append((read_id, distance, ref_id, read_seq, ref_seq, None))
    
    print(f"      Processing {len(all_work_items)} read-reference pairs using {threads} threads...")
    
    # Process in parallel
    position_errors = defaultdict(lambda: {'subs': 0, 'ins': 0, 'del': 0})
    
    with Pool(threads) as pool:
        # Process with progress bar
        results = list(tqdm(
            pool.imap(analyze_single_alignment, all_work_items, chunksize=128),
            total=len(all_work_items),
            desc=f"      Analyzing {data_type} alignments"
        ))
    
    # Aggregate results
    for error_data, distance, dist_for_analysis, errors_by_pos in results:
        if errors_by_pos:
            for pos, error_type in errors_by_pos:
                if error_type == 'substitution':
                    position_errors[pos]['subs'] += 1
                elif error_type == 'insertion':
                    position_errors[pos]['ins'] += 1
                elif error_type == 'deletion':
                    position_errors[pos]['del'] += 1
    
    print(f"      Found errors at {len(position_errors)} positions")
    return position_errors

def print_genome_counts(real_df, sim_df, dataset_name):
    """Print genome counts in TSV format"""
    # Count reads per genome
    real_counts = real_df.groupby('genome').size()
    sim_counts = sim_df.groupby('genome').size()
    
    # Get all genomes
    all_genomes = sorted(set(real_counts.index) | set(sim_counts.index))
    
    print(f"\n{dataset_name} Genome Counts:")
    print("Genome\tReal\tSimulated")
    for genome in all_genomes:
        real_count = real_counts.get(genome, 0)
        sim_count = sim_counts.get(genome, 0)
        print(f"{genome}\t{real_count}\t{sim_count}")



def generate_subread_accuracy_figure(output_dir, dpi=600):
    """Generate Subread Accuracy Optimization figure"""
    print("\nGenerating Subread Accuracy figure...")
    
    # Define paths to linesearch results - FIXED ORDER
    base_dirs = [
        ('Zymo/Titan', '/data/shoreline/Simulator_datasets/Zymo_Titan_D6300/linesearch_analysis'),
        ('ATCC/16S', '/data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/linesearch_analysis'),
        ('Phylotag/16S', '/data/shoreline/Simulator_datasets/Phylotag/linesearch_analysis')
    ]
    
    # Changed to 1 row, 3 columns
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    for idx, (dataset_name, linesearch_dir) in enumerate(base_dirs):
        ax = axes[idx]
        
        # Read the CSV file with all results
        csv_path = os.path.join(linesearch_dir, 'all_linesearch_results.csv')
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            df_sorted = df.sort_values('accuracy')
            
            # Plot KL divergence
            ax.plot(df_sorted['accuracy'], df_sorted['kl_symmetric'], 
                   'b-o', linewidth=2, markersize=6)
            
            # Highlight minimum
            min_idx = df_sorted['kl_symmetric'].idxmin()
            min_row = df_sorted.loc[min_idx]
            ax.plot(min_row['accuracy'], min_row['kl_symmetric'], 
                   'ro', markersize=10)
            
            # Add annotation for best value to ALL plots
            ax.annotate(f'Best: {min_row["accuracy"]:.2f}\nKL: {min_row["kl_symmetric"]:.4f}',
                       xy=(min_row['accuracy'], min_row['kl_symmetric']),
                       xytext=(10, 10), textcoords='offset points',
                       bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.7),
                       arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'),
                       fontweight='bold')
            
            # Customize plot
            ax.set_xlabel('Subread Accuracy', fontweight='bold')
            ax.set_ylabel('Symmetric KL Divergence' if idx == 0 else '', fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.set_xlim(0.45, 1.05)
            
            # Add dataset title above each subplot
            ax.set_title(dataset_name, fontsize=14, fontweight='bold')
            
            # Bold tick labels
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontweight('bold')
    
    plt.suptitle('Subread Accuracy Optimization', fontsize=18, fontweight='bold')
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, 'subread_accuracy_optimization.png')
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_path}")

def generate_errors_by_position_figure(output_dir, dpi=600, threads=None):
    """Generate Errors by Position figure"""
    print("\nGenerating Errors by Position figure...")
    
    # Define paths
    datasets = {
        'Zymo/Titan': '/data/shoreline/Simulator_datasets/Zymo_Titan_D6300/Analysis_final',
        'ATCC/16S': '/data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/Analysis_final',
        'Phylotag/16S': '/data/shoreline/Simulator_datasets/Phylotag/Analysis_final'
    }
    
    fig, axes = plt.subplots(3, 2, figsize=(16, 12))
    
    # Track the maximum position across all datasets for x-axis alignment
    max_positions = []
    
    # For printing error distribution tables
    error_tables = {}
    
    for idx, (dataset_name, base_dir) in enumerate(datasets.items()):
        error_tables[dataset_name] = {}
        
        for col_idx, data_type in enumerate(['real', 'simulated']):
            ax = axes[idx, col_idx]
            
            # Generate error position data on the fly
            try:
                if data_type == 'real':
                    closest_file = os.path.join(base_dir, 'ref-closest-real.txt')
                    reads_fasta = os.path.join(base_dir, 'demux/amplicons.fa')
                    reference_fasta = os.path.join(base_dir, 'ref-deduped.fa')
                else:  # simulated
                    closest_file = os.path.join(base_dir, 'simulated/ref-closest-sim.txt')
                    reads_fasta = os.path.join(base_dir, 'simulated/demux/amplicons.fa')
                    reference_fasta = os.path.join(base_dir, 'ref-deduped.fa')
                
                # Check if files exist
                if not all(os.path.exists(f) for f in [closest_file, reads_fasta, reference_fasta]):
                    print(f"  WARNING: Missing files for {dataset_name} {data_type}")
                    ax.text(0.5, 0.5, 'Data not available', ha='center', va='center', 
                           transform=ax.transAxes, fontsize=14, fontweight='bold', color='red')
                    ax.set_xlim(0, 1)
                    ax.set_ylim(0, 1)
                    continue
                
                # Generate position errors with parallel processing
                position_errors = generate_error_position_data_parallel(
                    closest_file, reads_fasta, reference_fasta, data_type, threads
                )
                
                if not position_errors:
                    print(f"  WARNING: No error data generated for {dataset_name} {data_type}")
                    ax.text(0.5, 0.5, 'No error data', ha='center', va='center', 
                           transform=ax.transAxes, fontsize=14, fontweight='bold', color='red')
                    ax.set_xlim(0, 1)
                    ax.set_ylim(0, 1)
                    continue
                
                # Process the data
                positions = sorted(position_errors.keys())
                subs_counts = [position_errors[pos]['subs'] for pos in positions]
                ins_counts = [position_errors[pos]['ins'] for pos in positions]
                del_counts = [position_errors[pos]['del'] for pos in positions]
                
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
                
                # Calculate error distribution by ranges
                max_pos = max(positions) if positions else 1600
                ranges = [(0, max_pos*0.1), (max_pos*0.1, max_pos*0.2), (max_pos*0.2, max_pos*0.3),
                         (max_pos*0.3, max_pos*0.4), (max_pos*0.4, max_pos*0.5), (max_pos*0.5, max_pos*0.6),
                         (max_pos*0.6, max_pos*0.7), (max_pos*0.7, max_pos*0.8), (max_pos*0.8, max_pos*0.9),
                         (max_pos*0.9, max_pos)]
                
                range_stats = []
                for r_start, r_end in ranges:
                    indices = [i for i, p in enumerate(positions) if r_start <= p < r_end]
                    if indices:
                        avg_subs = np.mean([subs_prop[i] for i in indices])
                        avg_ins = np.mean([ins_prop[i] for i in indices])
                        avg_del = np.mean([del_prop[i] for i in indices])
                        range_stats.append({
                            'range': f'{int(r_start)}-{int(r_end)}',
                            'subs': avg_subs,
                            'ins': avg_ins,
                            'del': avg_del
                        })
                
                error_tables[dataset_name][data_type] = range_stats
                
                # Track max position
                if positions:
                    max_positions.append(max(positions))
                
                # Stack plot - only add labels to first subplot
                if idx == 0 and col_idx == 0:
                    ax.stackplot(positions, subs_prop, ins_prop, del_prop,
                                labels=['Substitutions', 'Insertions', 'Deletions'],
                                colors=['#D62728', '#2CA02C', '#1F77B4'],
                                alpha=0.8)
                else:
                    ax.stackplot(positions, subs_prop, ins_prop, del_prop,
                                colors=['#D62728', '#2CA02C', '#1F77B4'],
                                alpha=0.8)
                
                ax.set_xlabel('Position in Alignment (bp)' if idx == 2 else '', fontweight='bold')
                ax.set_ylabel('Proportion of Error Type' if col_idx == 0 else '', fontweight='bold')
                ax.set_ylim(0, 1)
                
                # Set x-axis to fit data
                if positions:
                    ax.set_xlim(0, max(positions))
                
                ax.grid(True, alpha=0.3)
                
                # Add data type label
                data_label = 'Real' if col_idx == 0 else 'Simulated'
                ax.text(0.98, 0.98, data_label, transform=ax.transAxes,
                       fontsize=12, va='top', ha='right', fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
                
                # Bold tick labels
                for label in ax.get_xticklabels() + ax.get_yticklabels():
                    label.set_fontweight('bold')
                
                # Only show legend once
                if idx == 0 and col_idx == 0:
                    ax.legend(loc='upper left', prop={'weight': 'bold'})
                    
            except Exception as e:
                print(f"  ERROR processing {dataset_name} {data_type}: {str(e)}")
                ax.text(0.5, 0.5, 'Processing error', ha='center', va='center', 
                       transform=ax.transAxes, fontsize=14, fontweight='bold', color='red')
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
        
        # Add dataset label
        fig.text(0.01, 0.83 - idx * 0.33, dataset_name, 
                fontsize=14, fontweight='bold', rotation=90, va='center')
    
    plt.suptitle('Error Composition by Position', fontsize=18, fontweight='bold')
    plt.tight_layout()
    plt.subplots_adjust(left=0.08)
    
    output_path = os.path.join(output_dir, 'errors_by_position.png')
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_path}")
    
    # Print error distribution tables
    print("\nError Distribution by Position Range:")
    for dataset_name, data_dict in error_tables.items():
        for data_type, stats in data_dict.items():
            if stats:
                print(f"\n{dataset_name} - {data_type}:")
                print("Range\tSubs%\tIns%\tDel%")
                for s in stats:
                    print(f"{s['range']}\t{s['subs']*100:.1f}\t{s['ins']*100:.1f}\t{s['del']*100:.1f}")

def calculate_expected_abundances_corrected(dedup_mapping_file, genome_abunds_file):
    """
    Calculate expected abundances correctly accounting for ASV copy numbers.
    
    The genome abundance TSV gives molarity-based proportions.
    After amplification, each ASV contributes equally to the pool.
    So we need to weight by the number of ASVs per genome.
    """
    print("  Calculating corrected expected abundances...")
    
    # Read genome abundances (molarity-based)
    genome_abunds = pd.read_csv(genome_abunds_file, sep='\t')
    molarity_dict = {}
    for _, row in genome_abunds.iterrows():
        genome = row['GenomeID']
        molarity = row['Proportion'] / 100.0  # Convert percentage to proportion
        molarity_dict[genome] = molarity
    
    # Read dedup mapping and count ACTUAL ASVs (using num_original)
    dedup_df = pd.read_csv(dedup_mapping_file, sep='\t')
    
    # Count total ASVs per genome (sum of num_original)
    asvs_per_genome = dedup_df.groupby('genome')['num_original'].sum().to_dict()
    
    # Get ASV details per genome
    asvs_by_genome = {}
    for genome in dedup_df['genome'].unique():
        genome_df = dedup_df[dedup_df['genome'] == genome]
        asvs_by_genome[genome] = [(row['dedup_id'], row['num_original']) 
                                  for _, row in genome_df.iterrows()]
    
    # Print detailed ASV information
    print("\n  ASVs per genome (counting all copies):")
    print("  " + "-" * 60)
    total_asvs = 0
    for genome in sorted(asvs_per_genome.keys()):
        count = asvs_per_genome[genome]
        total_asvs += count
        print(f"  {genome}: {count} ASVs")
        # Print each deduplicated ASV ID with its copy number
        for asv_id, num_copies in asvs_by_genome[genome]:
            if num_copies > 1:
                print(f"    - {asv_id} (represents {num_copies} identical copies)")
            else:
                print(f"    - {asv_id}")
    print("  " + "-" * 60)
    print(f"  TOTAL ASVs: {total_asvs}")
    print()
    
    # Calculate expected abundance after amplification
    # Each ASV gets equal representation, so we weight by ASV count
    weighted_abundances = {}
    total_weighted = 0
    
    for genome, n_asvs in asvs_per_genome.items():
        if genome in molarity_dict:
            # Weighted abundance = molarity * number of ASVs
            weighted = molarity_dict[genome] * n_asvs
            weighted_abundances[genome] = weighted
            total_weighted += weighted
    
    # Normalize to sum to 1
    expected_abundances = {}
    for genome, weighted in weighted_abundances.items():
        expected_abundances[genome] = weighted / total_weighted if total_weighted > 0 else 0
    
    # Print summary
    print(f"    Genomes with molarity data: {len(molarity_dict)}")
    print(f"    Genomes with ASV counts: {len(asvs_per_genome)}")
    print(f"    Total original ASVs across all genomes: {total_asvs}")
    
    return expected_abundances, asvs_per_genome

def generate_figure1_violins(output_dir, dpi=600):
    """Generate Figure 1: Multi-panel violin plots"""
    print("Generating Figure 1: Multi-panel violin plots...")
    
    # Define paths - FIXED ORDER
    datasets = [
        ('Zymo/Titan', {
            'real': "/data/shoreline/Simulator_datasets/Zymo_Titan_D6300/Analysis_final/ref-closest-real.txt",
            'sim': "/data/shoreline/Simulator_datasets/Zymo_Titan_D6300/Analysis_final/simulated/ref-closest-sim.txt"
        }),
        ('ATCC/16S', {
            'real': "/data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/Analysis_final/ref-closest-real.txt",
            'sim': "/data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/Analysis_final/simulated/ref-closest-sim.txt"
        }),
        ('Phylotag/16S', {
            'real': "/data/shoreline/Simulator_datasets/Phylotag/Analysis_final/ref-closest-real.txt",
            'sim': "/data/shoreline/Simulator_datasets/Phylotag/Analysis_final/simulated/ref-closest-sim.txt"
        })
    ]
    
    # Read data
    all_data = []
    for dataset_name, paths in datasets:
        real_df = read_closest(paths['real'])
        sim_df = read_closest(paths['sim'])
        all_data.append((dataset_name, real_df, sim_df))
    
    # Create figure
    fig = plt.figure(figsize=(16, 20))
    gs = GridSpec(3, 1, figure=fig, height_ratios=[1, 1, 1], hspace=0.35)
    
    axes = [fig.add_subplot(gs[i, 0]) for i in range(3)]
    
    # Process each dataset
    for idx, (dataset_name, real_df, sim_df) in enumerate(all_data):
        ax = axes[idx]
        
        # Add genome column
        real_df = real_df.copy()
        sim_df = sim_df.copy()
        real_df['genome'] = real_df['reference_id'].apply(extract_genome)
        sim_df['genome'] = sim_df['reference_id'].apply(extract_genome)
        
        # Print genome counts for this dataset
        print_genome_counts(real_df, sim_df, dataset_name)
        
        # Format genome names with italics
        real_df['genome_display'] = real_df['genome'].apply(format_genome_name)
        sim_df['genome_display'] = sim_df['genome'].apply(format_genome_name)
        
        # Add source column
        real_df['source'] = 'Real'
        sim_df['source'] = 'Simulated'
        
        # Combine datasets
        combined = pd.concat([real_df, sim_df])
        
        # Get all unique genomes and sort
        all_genomes_display = sorted(combined['genome_display'].unique())
        genome_data = combined[combined['genome_display'].isin(all_genomes_display)]
        
        if len(genome_data) > 0:
            # Create violin plot
            violin_parts = sns.violinplot(
                x='genome_display', 
                y='edit_distance', 
                hue='source', 
                data=genome_data, 
                split=True, 
                inner='quart', 
                palette=['blue', 'red'],
                ax=ax,
                order=all_genomes_display
            )
            
            # Customize the plot
            ax.set_xlabel('Genome' if idx == 2 else '', fontweight='bold')
            ax.set_ylabel('Edit Distance', fontweight='bold')
            
            # Rotate labels with fixed ticklabels
            labels = [label.get_text() for label in ax.get_xticklabels()]
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, rotation=45, ha='right', fontweight='bold')
            
            # Only show legend on second panel
            if idx == 1:
                ax.legend(title='Data Source', loc='upper right', fontsize=10, 
                         title_fontsize=11, prop={'weight': 'bold'})
            else:
                ax.get_legend().remove()
            
            # Add grid
            ax.grid(True, axis='y', alpha=0.3)
            ax.set_axisbelow(True)
            
            # Add dataset label - CONSISTENT POSITIONING for all datasets
            x_position = -0.07  # Fixed position for all datasets
            ax.text(x_position, 0.5, dataset_name, transform=ax.transAxes,
                    fontsize=16, fontweight='bold', va='center', ha='center',
                    rotation=90)
            
            # Add sample counts
            real_count = len(real_df)
            sim_count = len(sim_df)
            count_text = f'n={real_count:,} real, {sim_count:,} sim'
            ax.text(0.02, 0.98, count_text, transform=ax.transAxes,
                    fontsize=9, va='top', ha='left', fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
            
            # Bold tick labels
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontweight('bold')
    
    # Add main title
    fig.suptitle('Edit Distance Distributions by Genome', 
                 fontsize=20, fontweight='bold', y=0.98)
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(top=0.96, left=0.1, right=0.98, bottom=0.05)
    
    # Save figure
    output_path = os.path.join(output_dir, 'multi_panel_violins.png')
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_path}")
    
    # Print note about violin plots
    print("\nNote: Violin plots show probability density distributions, not raw counts.")
    print("The area under each violin represents the distribution shape, normalized for visualization.")
    print("The split violins allow visual comparison of distributions between real and simulated data.")


def generate_abundance_figure(output_dir, dpi=600):
    """Generate Abundance Estimation figure"""
    print("\nGenerating Abundance Estimation figure...")
    
    # Define paths - FIXED ORDER
    datasets = [
        ('Zymo/Titan', {
            'plot_dir': '/data/shoreline/Simulator_datasets/Zymo_Titan_D6300/Analysis_final/abundance_plots_genome',
            'dedup_mapping': '/data/shoreline/Simulator_datasets/Zymo_Titan_D6300/Analysis_final/ref-dedup-mapping.tsv',
            'genome_abunds': '/data/shoreline/Simulator_datasets/Zymo_Titan_D6300/sim/amplicons/genome_abunds.tsv'
        }),
        ('ATCC/16S', {
            'plot_dir': '/data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/Analysis_final/abundance_plots_genome',
            'dedup_mapping': '/data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/Analysis_final/ref-dedup-mapping.tsv',
            'genome_abunds': '/data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/genomes/genome_abunds.tsv'
        }),
        ('Phylotag/16S', {
            'plot_dir': '/data/shoreline/Simulator_datasets/Phylotag/Analysis_final/abundance_plots_genome',
            'dedup_mapping': '/data/shoreline/Simulator_datasets/Phylotag/Analysis_final/ref-dedup-mapping.tsv',
            'genome_abunds': '/data/shoreline/Simulator_datasets/Phylotag/genome_abunds.tsv'
        })
    ]
    
    # Changed to 3x1 layout (3 rows, 1 column)
    fig, axes = plt.subplots(3, 1, figsize=(12, 20))
    
    for idx, (dataset_name, paths) in enumerate(datasets):
        ax = axes[idx]
        
        # Calculate corrected expected abundances
        expected_corrected, asvs_per_genome = calculate_expected_abundances_corrected(
            paths['dedup_mapping'], paths['genome_abunds']
        )
        
        # Load the observed data CSV
        csv_path = os.path.join(paths['plot_dir'], 'expected_vs_observed_data_genome.csv')
        if os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
            
            # Replace expected values with corrected ones
            df['Expected_Original'] = df['Expected']  # Keep original for comparison
            df['Expected'] = df['ID'].map(expected_corrected).fillna(0)
            
            # Sort by corrected expected abundance
            df = df.sort_values('Expected', ascending=False)
            
            # Create bar plot
            x = np.arange(len(df))
            width = 0.25
            
            # Plot bars - add labels to all bars for legend purposes
            bars1 = ax.bar(x - width, df['Expected'], width, 
                           label='Expected', color='blue', alpha=0.8)
            bars2 = ax.bar(x, df['Real Observed'], width, 
                           label='Real Observed', color='green', alpha=0.8)
            bars3 = ax.bar(x + width, df['Simulated Observed'], width, 
                           label='Simulated Observed', color='red', alpha=0.8)
            
            # Customize plot
            ax.set_xlabel('Genome' if idx == 2 else '', fontweight='bold')
            ax.set_ylabel('Relative Abundance', fontweight='bold')
            ax.set_title(dataset_name, fontsize=14, fontweight='bold')
            ax.set_xticks(x)
            
            # Format x-axis labels with italics
            labels = [format_genome_name(genome) for genome in df['ID']]
            ax.set_xticklabels(labels, rotation=45, ha='right', fontweight='bold')
            ax.grid(axis='y', alpha=0.3)
            
            # Apply log scale for ATCC and Phylotag only
            if dataset_name in ['ATCC/16S', 'Phylotag/16S']:
                ax.set_yscale('log')
                # Set reasonable y-axis limits for log scale
                y_min = min(df[['Expected', 'Real Observed', 'Simulated Observed']].min().min(), 1e-4)
                y_max = df[['Expected', 'Real Observed', 'Simulated Observed']].max().max() * 2
                ax.set_ylim(y_min, y_max)
            
            # Calculate correlations with corrected expected values
            real_corr = np.corrcoef(df['Expected'], df['Real Observed'])[0, 1]
            sim_corr = np.corrcoef(df['Expected'], df['Simulated Observed'])[0, 1]
            
            # Add correlation text
            corr_text = f'Real R²={real_corr**2:.3f}\nSim R²={sim_corr**2:.3f}'
            ax.text(0.98, 0.98, corr_text, transform=ax.transAxes,
                   fontsize=10, va='top', ha='right', fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
            
            # Bold tick labels
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontweight('bold')
            
            # Show legend on the second panel
            if idx == 1:
                ax.legend(loc='center right', prop={'weight': 'bold'})
    
    plt.suptitle('Expected vs Observed Abundances', fontsize=18, fontweight='bold')
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, 'abundance_estimation.png')
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_path}")

def generate_length_distribution_figure(output_dir, dpi=600):
    """Generate Length Distribution figure"""
    print("\nGenerating Length Distribution figure...")
    
    # Define paths to read FASTQ files and x-axis ranges - FIXED ORDER
    datasets = [
        ('Zymo/Titan', {
            'real': '/data/shoreline/Simulator_datasets/Zymo_Titan_D6300/real/m54215_200618_132850.Q20.fastq',
            'sim': '/data/shoreline/Simulator_datasets/Zymo_Titan_D6300/Analysis_final/simulated/combined_reads.fastq',
            'xlim': (2000, 3000)
        }),
        ('ATCC/16S', {
            'real': '/data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/All_reads.fastq',
            'sim': '/data/shoreline/Simulator_datasets/ATCC_16S_MSA-1003/Analysis_final/simulated/combined_reads.fastq',
            'xlim': (1000, 2000)
        }),
        ('Phylotag/16S', {
            'real': '/data/shoreline/Simulator_datasets/Phylotag/real/Mock1-5.fastq',
            'sim': '/data/shoreline/Simulator_datasets/Phylotag/Analysis_final/simulated/combined_reads.fastq',
            'xlim': (1000, 2000)
        })
    ]
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    for idx, (dataset_name, paths) in enumerate(datasets):
        ax = axes[idx]
        
        # Read lengths from FASTQ files
        real_lengths = []
        sim_lengths = []
        
        # Read real data
        if os.path.exists(paths['real']):
            print(f"  Reading real lengths from {paths['real']}")
            with open(paths['real'], 'r') as f:
                line_count = 0
                for line in f:
                    if line_count % 4 == 1:  # Sequence line
                        real_lengths.append(len(line.strip()))
                    line_count += 1
                    if len(real_lengths) >= 10000:  # Sample for speed
                        break
        
        # Read simulated data
        if os.path.exists(paths['sim']):
            print(f"  Reading simulated lengths from {paths['sim']}")
            with open(paths['sim'], 'r') as f:
                line_count = 0
                for line in f:
                    if line_count % 4 == 1:  # Sequence line
                        sim_lengths.append(len(line.strip()))
                    line_count += 1
                    if len(sim_lengths) >= 10000:  # Sample for speed
                        break
        
        # Create overlaid histograms
        if real_lengths and sim_lengths:
            # Use the fixed x-axis range for bins
            xmin, xmax = paths['xlim']
            bins = np.linspace(xmin, xmax, 50)
            
            # Plot histograms
            ax.hist(real_lengths, bins=bins, alpha=0.5, label='Real' if idx == 0 else None, 
                    color='blue', density=True)
            ax.hist(sim_lengths, bins=bins, alpha=0.5, label='Simulated' if idx == 0 else None, 
                    color='red', density=True)
            
            # Set x-axis limits
            ax.set_xlim(xmin, xmax)
        
        # Customize plot
        ax.set_xlabel('Read Length (bp)', fontweight='bold')
        ax.set_ylabel('Density' if idx == 0 else '', fontweight='bold')
        ax.set_title(dataset_name, fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # Bold tick labels
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')
        
        # Only show legend on first panel
        if idx == 0:
            ax.legend(loc='upper right', prop={'weight': 'bold'})
    
    plt.suptitle('Read Length Distributions', fontsize=18, fontweight='bold')
    plt.tight_layout()
    
    output_path = os.path.join(output_dir, 'Figure_length_distributions.png')
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_path}")

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print(f"Generating all publication figures...")
    print(f"Output directory: {args.output_dir}")
    print(f"DPI: {args.dpi}")
    print(f"Threads: {args.threads if args.threads else 'all available CPUs'}")
    print("")
    
    # Generate all figures
    generate_figure1_violins(args.output_dir, args.dpi)
    generate_subread_accuracy_figure(args.output_dir, args.dpi)
    generate_errors_by_position_figure(args.output_dir, args.dpi, args.threads)
    generate_abundance_figure(args.output_dir, args.dpi)
    generate_length_distribution_figure(args.output_dir, args.dpi)
    
    print("\nAll figures generated successfully!")

if __name__ == '__main__':
    main()
