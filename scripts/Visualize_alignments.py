#!/usr/bin/env python3
"""
Visualize_alignments.py
Create visual representations of multiple sequence alignments from FASTA files.
Generates both text-based and graphical visualizations.
"""

import argparse
import random
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from collections import Counter
import seaborn as sns
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(description='Visualize multiple sequence alignments')
    parser.add_argument('-i', '--input_dir', required=True,
                        help='Directory containing .aligned files')
    parser.add_argument('-o', '--output_dir', required=True,
                        help='Output directory for visualizations')
    parser.add_argument('-f', '--format', choices=['png', 'pdf', 'svg'], default='png',
                        help='Output format for plots (default: png)')
    parser.add_argument('--max_seqs', type=int, default=50,
                        help='Maximum sequences to show in detailed view (default: 50)')
    parser.add_argument('--window_size', type=int, default=100,
                        help='Window size for conservation plot (default: 100)')
    parser.add_argument('--text_output', action='store_true',
                        help='Also generate text-based alignment views')
    return parser.parse_args()

def read_alignment(filepath):
    """Read alignment from FASTA file."""
    try:
        # Read as FASTA alignment
        alignment = AlignIO.read(filepath, "fasta")
        return alignment
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None

def calculate_conservation(alignment, window_size=1):
    """Calculate conservation score for each position."""
    conservation = []
    align_array = np.array([list(str(rec.seq)) for rec in alignment])
    
    for i in range(0, alignment.get_alignment_length(), window_size):
        window = align_array[:, i:i+window_size]
        # Calculate conservation as 1 - Shannon entropy
        scores = []
        for j in range(window.shape[1]):
            col = window[:, j]
            # Skip gap-only columns
            non_gaps = col[col != '-']
            if len(non_gaps) == 0:
                scores.append(0)
                continue
            
            # Calculate entropy
            counts = Counter(non_gaps)
            total = sum(counts.values())
            entropy = -sum((c/total) * np.log2(c/total) for c in counts.values() if c > 0)
            max_entropy = np.log2(min(4, len(set(non_gaps))))  # Max entropy for DNA
            
            # Normalize to 0-1 scale
            conservation_score = 1 - (entropy / max_entropy) if max_entropy > 0 else 1
            scores.append(conservation_score)
        
        # Average conservation in window
        conservation.extend([np.mean(scores)] * min(window_size, window.shape[1]))
    
    return conservation[:alignment.get_alignment_length()]

def create_alignment_heatmap(alignment, output_path, title="", max_seqs=50):
    """Create a heatmap visualization of the alignment."""
    # Limit number of sequences for readability
    if len(alignment) > max_seqs:
        total_seqs=len(alignment)
        alignment = MultipleSeqAlignment(random.sample(alignment,max_seqs))
        title += f" (showing {max_seqs}/{total_seqs} sequences)"
    
    # Convert to numeric matrix
    base_to_num = {'A': 1, 'T': 2, 'G': 3, 'C': 4, '-': 0, 'N': 5}
    
    matrix = []
    for record in alignment:
        seq_nums = [base_to_num.get(base.upper(), 5) for base in str(record.seq)]
        matrix.append(seq_nums)
    
    matrix = np.array(matrix)
    
    # Create figure
    fig_height = max(6, len(alignment) * 0.3)
    fig_width = max(12, alignment.get_alignment_length() * 0.05)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(fig_width, fig_height + 2),
                                    gridspec_kw={'height_ratios': [1, 5]})
    
    # Conservation plot
    conservation = calculate_conservation(alignment)
    ax1.plot(conservation, color='darkblue', linewidth=1)
    ax1.fill_between(range(len(conservation)), conservation, alpha=0.3, color='darkblue')
    ax1.set_xlim(0, len(conservation))
    ax1.set_ylim(0, 1)
    ax1.set_ylabel('Conservation')
    ax1.set_title(f'Alignment Conservation - {title}')
    ax1.grid(True, alpha=0.3)
    
    # Alignment heatmap
    cmap = plt.cm.colors.ListedColormap(['white', 'red', 'blue', 'green', 'yellow', 'gray'])
    im = ax2.imshow(matrix, cmap=cmap, aspect='auto', interpolation='nearest')
    
    # Add sequence names
    ax2.set_yticks(range(len(alignment)))
    ax2.set_yticklabels([rec.id[:30] + "..." if len(rec.id) > 30 else rec.id 
                         for rec in alignment], fontsize=8)
    
    # Add position labels
    positions = list(range(0, alignment.get_alignment_length(), 100))
    ax2.set_xticks(positions)
    ax2.set_xticklabels(positions, rotation=45)
    ax2.set_xlabel('Position')
    
    # Add colorbar legend
    cbar = plt.colorbar(im, ax=ax2, orientation='horizontal', pad=0.05, shrink=0.5)
    cbar.set_ticks([0, 1, 2, 3, 4, 5])
    cbar.set_ticklabels(['-', 'A', 'T', 'G', 'C', 'N'])
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

def create_gap_distribution_plot(alignment, output_path, title=""):
    """Create a plot showing gap distribution across the alignment."""
    gap_counts = []
    positions = []
    
    for i in range(alignment.get_alignment_length()):
        col = alignment[:, i]
        gap_count = col.count('-')
        gap_counts.append(gap_count)
        positions.append(i)
    
    plt.figure(figsize=(12, 6))
    
    # Gap count plot
    plt.subplot(2, 1, 1)
    plt.plot(positions, gap_counts, color='darkred', linewidth=1)
    plt.fill_between(positions, gap_counts, alpha=0.3, color='darkred')
    plt.xlabel('Position')
    plt.ylabel('Number of Gaps')
    plt.title(f'Gap Distribution - {title}')
    plt.grid(True, alpha=0.3)
    
    # Gap percentage plot
    plt.subplot(2, 1, 2)
    gap_percentages = [count / len(alignment) * 100 for count in gap_counts]
    plt.plot(positions, gap_percentages, color='darkblue', linewidth=1)
    plt.fill_between(positions, gap_percentages, alpha=0.3, color='darkblue')
    plt.xlabel('Position')
    plt.ylabel('Gap Percentage (%)')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

def create_text_alignment_view(alignment, output_path, max_seqs=20, max_width=80):
    """Create a text-based view of the alignment."""
    with open(output_path, 'w') as f:
        # Limit sequences
        display_seqs = min(len(alignment), max_seqs)
        
        # Write in blocks
        for start in range(0, alignment.get_alignment_length(), max_width):
            end = min(start + max_width, alignment.get_alignment_length())
            
            # Write position ruler
            f.write(" " * 30 + "".join(f"{i%10}" for i in range(start, end)) + "\n")
            
            # Write sequences
            for i in range(display_seqs):
                seq_id = alignment[i].id[:25].ljust(25)
                seq_slice = str(alignment[i].seq[start:end])
                f.write(f"{seq_id} {seq_slice}\n")
            
            if display_seqs < len(alignment):
                f.write(f"... and {len(alignment) - display_seqs} more sequences ...\n")
            
            # Write consensus line
            consensus = []
            for pos in range(start, end):
                col = alignment[:, pos]
                non_gaps = [b for b in col if b != '-']
                if not non_gaps:
                    consensus.append('-')
                else:
                    counter = Counter(non_gaps)
                    most_common = counter.most_common(1)[0]
                    if most_common[1] / len(non_gaps) > 0.7:
                        consensus.append(most_common[0].upper())
                    elif most_common[1] / len(non_gaps) > 0.5:
                        consensus.append(most_common[0].lower())
                    else:
                        consensus.append('.')
            
            f.write("Consensus".ljust(25) + " " + "".join(consensus) + "\n")
            f.write("\n")

def create_summary_plot(all_stats, output_path):
    """Create a summary plot of all alignments."""
    if not all_stats:
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Number of sequences histogram
    ax = axes[0, 0]
    num_seqs = [s['num_seqs'] for s in all_stats]
    ax.hist(num_seqs, bins=20, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Number of Sequences')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Alignment Sizes')
    ax.grid(True, alpha=0.3)
    
    # Alignment length histogram
    ax = axes[0, 1]
    lengths = [s['length'] for s in all_stats]
    ax.hist(lengths, bins=20, edgecolor='black', alpha=0.7, color='green')
    ax.set_xlabel('Alignment Length (bp)')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Alignment Lengths')
    ax.grid(True, alpha=0.3)
    
    # Average conservation
    ax = axes[1, 0]
    conservations = [s['avg_conservation'] for s in all_stats]
    ax.hist(conservations, bins=20, edgecolor='black', alpha=0.7, color='orange')
    ax.set_xlabel('Average Conservation Score')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Conservation Scores')
    ax.grid(True, alpha=0.3)
    
    # Gap percentage
    ax = axes[1, 1]
    gap_percentages = [s['gap_percentage'] for s in all_stats]
    ax.hist(gap_percentages, bins=20, edgecolor='black', alpha=0.7, color='red')
    ax.set_xlabel('Gap Percentage (%)')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Gap Content')
    ax.grid(True, alpha=0.3)
    
    plt.suptitle('Alignment Statistics Summary', fontsize=16)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

def process_alignment_file(filepath, output_dir, args):
    """Process a single alignment file."""
    basename = os.path.basename(filepath).replace('.aligned', '')
    
    # Read alignment
    alignment = read_alignment(filepath)
    if alignment is None:
        return None
    
    # Calculate statistics
    stats = {
        'name': basename,
        'num_seqs': len(alignment),
        'length': alignment.get_alignment_length(),
        'avg_conservation': np.mean(calculate_conservation(alignment)),
        'gap_percentage': sum(str(rec.seq).count('-') for rec in alignment) / 
                          (len(alignment) * alignment.get_alignment_length()) * 100
    }
    
    # Create visualizations
    # 1. Alignment heatmap
    heatmap_path = os.path.join(output_dir, f"{basename}_heatmap.{args.format}")
    create_alignment_heatmap(alignment, heatmap_path, title=basename, max_seqs=args.max_seqs)
    
    # 2. Gap distribution
    gap_path = os.path.join(output_dir, f"{basename}_gaps.{args.format}")
    create_gap_distribution_plot(alignment, gap_path, title=basename)
    
    # 3. Text view if requested
    if args.text_output:
        text_path = os.path.join(output_dir, f"{basename}_alignment.txt")
        create_text_alignment_view(alignment, text_path)
    
    return stats

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Find all alignment files
    alignment_files = list(Path(args.input_dir).glob('*.aligned'))
    print(f"Found {len(alignment_files)} alignment files")
    
    # Process each alignment
    all_stats = []
    for i, filepath in enumerate(alignment_files, 1):
        print(f"Processing {i}/{len(alignment_files)}: {filepath.name}")
        stats = process_alignment_file(str(filepath), args.output_dir, args)
        if stats:
            all_stats.append(stats)
    
    # Create summary plot
    if all_stats:
        summary_path = os.path.join(args.output_dir, f"alignment_summary.{args.format}")
        create_summary_plot(all_stats, summary_path)
        
        # Save statistics to CSV
        import pandas as pd
        df = pd.DataFrame(all_stats)
        df.to_csv(os.path.join(args.output_dir, "alignment_statistics.csv"), index=False)
        print(f"\nSaved statistics for {len(all_stats)} alignments")
    
    print(f"Visualizations saved to {args.output_dir}")

if __name__ == '__main__':
    main()
