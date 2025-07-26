#!/usr/bin/env python3
"""
Simple PacBio CCS FASTQ Comparison Script
Compares simulated vs real PacBio CCS data using only FASTQ information
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
import subprocess
import re
from pathlib import Path

def run_command(cmd):
    """Execute shell command and return output"""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout, result.stderr

def parse_fastq(fastq_file, max_reads=None):
    """Parse FASTQ file and extract basic metrics"""
    lengths = []
    qualities = []
    gc_contents = []
    sequences = []
    
    count = 0
    with open(fastq_file, 'r') as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            
            # Calculate metrics
            lengths.append(len(seq))
            qualities.append(np.mean([ord(q) - 33 for q in qual]))
            
            gc_count = seq.upper().count('G') + seq.upper().count('C')
            gc_contents.append((gc_count / len(seq)) * 100 if len(seq) > 0 else 0)
            
            sequences.append(seq)
            
            count += 1
            if max_reads and count >= max_reads:
                break
    
    return {
        'lengths': lengths,
        'qualities': qualities,
        'gc_contents': gc_contents,
        'sequences': sequences
    }

def calculate_n50(lengths):
    """Calculate N50 from list of lengths"""
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    cumsum = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= total_length / 2:
            return length
    return 0

def analyze_homopolymers(sequences, max_length=15):
    """Analyze homopolymer distributions"""
    homo_counts = defaultdict(int)
    
    for seq in sequences:
        for base in ['A', 'T', 'G', 'C']:
            pattern = f"{base}+"
            matches = re.finditer(pattern, seq.upper())
            for match in matches:
                length = len(match.group())
                if length >= 2 and length <= max_length:
                    homo_counts[length] += 1
    
    return homo_counts

def plot_comparison(real_data, sim_data, output_dir):
    """Create all comparison plots"""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Set style
    plt.style.use('seaborn-v0_8-darkgrid')
    
    # 1. Length distributions
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Violin plot
    data_list = []
    for length in real_data['lengths']:
        data_list.append({'Dataset': 'Real', 'Read Length': length})
    for length in sim_data['lengths']:
        data_list.append({'Dataset': 'Simulated', 'Read Length': length})
    
    df = pd.DataFrame(data_list)
    sns.violinplot(data=df, x='Dataset', y='Read Length', ax=ax1)
    ax1.set_title('Read Length Distribution')
    ax1.set_ylabel('Read Length (bp)')
    
    # Cumulative length plot
    for label, data, color in [('Real', real_data, 'blue'), ('Simulated', sim_data, 'orange')]:
        sorted_lengths = sorted(data['lengths'], reverse=True)
        cumsum = np.cumsum(sorted_lengths)
        ax2.plot(range(len(sorted_lengths)), cumsum, label=label, alpha=0.8, color=color)
    
    ax2.set_xlabel('Read Index')
    ax2.set_ylabel('Cumulative Length (bp)')
    ax2.set_title('Cumulative Read Length')
    ax2.legend()
    ax2.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'length_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Quality distributions
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Quality histogram
    ax1.hist(real_data['qualities'], bins=50, alpha=0.6, label='Real', density=True, color='blue')
    ax1.hist(sim_data['qualities'], bins=50, alpha=0.6, label='Simulated', density=True, color='orange')
    ax1.set_xlabel('Mean Read Quality Score')
    ax1.set_ylabel('Density')
    ax1.set_title('Read Quality Distribution')
    ax1.legend()
    
    # Quality vs Length scatter (subsample for performance)
    n_points = min(5000, len(real_data['lengths']))
    
    # Real data
    indices = np.random.choice(len(real_data['lengths']), n_points, replace=False)
    ax2.scatter([real_data['lengths'][i] for i in indices],
                [real_data['qualities'][i] for i in indices],
                alpha=0.3, s=1, label='Real', color='blue')
    
    # Simulated data
    indices = np.random.choice(len(sim_data['lengths']), n_points, replace=False)
    ax2.scatter([sim_data['lengths'][i] for i in indices],
                [sim_data['qualities'][i] for i in indices],
                alpha=0.3, s=1, label='Simulated', color='orange')
    
    ax2.set_xlabel('Read Length (bp)')
    ax2.set_ylabel('Mean Quality Score')
    ax2.set_title('Quality vs Length')
    ax2.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'quality_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. GC content distribution
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    
    ax.hist(real_data['gc_contents'], bins=50, alpha=0.6, label='Real', density=True, color='blue')
    ax.hist(sim_data['gc_contents'], bins=50, alpha=0.6, label='Simulated', density=True, color='orange')
    ax.set_xlabel('GC Content (%)')
    ax.set_ylabel('Density')
    ax.set_title('GC Content Distribution')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(output_dir / 'gc_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 4. Homopolymer distributions
    real_homo = analyze_homopolymers(real_data['sequences'])
    sim_homo = analyze_homopolymers(sim_data['sequences'])
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    max_length = max(max(real_homo.keys(), default=2), max(sim_homo.keys(), default=2))
    lengths = list(range(2, min(max_length + 1, 16)))
    
    real_counts = [real_homo.get(l, 0) for l in lengths]
    sim_counts = [sim_homo.get(l, 0) for l in lengths]
    
    x = np.arange(len(lengths))
    width = 0.35
    
    ax.bar(x - width/2, real_counts, width, label='Real', alpha=0.8, color='blue')
    ax.bar(x + width/2, sim_counts, width, label='Simulated', alpha=0.8, color='orange')
    
    ax.set_xlabel('Homopolymer Length')
    ax.set_ylabel('Count')
    ax.set_title('Homopolymer Length Distribution')
    ax.set_xticks(x)
    ax.set_xticklabels(lengths)
    ax.legend()
    ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'homopolymer_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()

def save_statistics(real_data, sim_data, output_dir):
    """Save comparison statistics to TSV files"""
    output_dir = Path(output_dir)
    
    # Basic statistics
    stats = {
        'Metric': ['Number of reads', 'Total bases', 'Mean length', 'Median length', 
                   'Min length', 'Max length', 'N50', 'Mean quality', 'Mean GC%'],
        'Real': [
            len(real_data['lengths']),
            sum(real_data['lengths']),
            np.mean(real_data['lengths']),
            np.median(real_data['lengths']),
            min(real_data['lengths']),
            max(real_data['lengths']),
            calculate_n50(real_data['lengths']),
            np.mean(real_data['qualities']),
            np.mean(real_data['gc_contents'])
        ],
        'Simulated': [
            len(sim_data['lengths']),
            sum(sim_data['lengths']),
            np.mean(sim_data['lengths']),
            np.median(sim_data['lengths']),
            min(sim_data['lengths']),
            max(sim_data['lengths']),
            calculate_n50(sim_data['lengths']),
            np.mean(sim_data['qualities']),
            np.mean(sim_data['gc_contents'])
        ]
    }
    
    stats_df = pd.DataFrame(stats)
    stats_df.to_csv(output_dir / 'basic_statistics.tsv', sep='\t', index=False)
    
    # Per-read metrics (first 10000 reads)
    n_reads = min(10000, len(real_data['lengths']), len(sim_data['lengths']))
    per_read_data = []
    
    for i in range(n_reads):
        if i < len(real_data['lengths']):
            per_read_data.append({
                'Dataset': 'Real',
                'Read_index': i,
                'Length': real_data['lengths'][i],
                'Mean_quality': real_data['qualities'][i],
                'GC_content': real_data['gc_contents'][i]
            })
        
        if i < len(sim_data['lengths']):
            per_read_data.append({
                'Dataset': 'Simulated',
                'Read_index': i,
                'Length': sim_data['lengths'][i],
                'Mean_quality': sim_data['qualities'][i],
                'GC_content': sim_data['gc_contents'][i]
            })
    
    per_read_df = pd.DataFrame(per_read_data)
    per_read_df.to_csv(output_dir / 'per_read_metrics.tsv', sep='\t', index=False)

def generate_report(output_dir):
    """Generate simple HTML report"""
    output_dir = Path(output_dir)
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>PacBio CCS FASTQ Comparison Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; max-width: 1200px; margin: auto; }}
            h1 {{ color: #333; text-align: center; }}
            h2 {{ color: #666; margin-top: 30px; }}
            img {{ max-width: 100%; height: auto; margin: 10px 0; display: block; }}
            .section {{ margin-bottom: 30px; }}
            table {{ border-collapse: collapse; width: 100%; margin-top: 10px; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
        </style>
    </head>
    <body>
        <h1>PacBio CCS FASTQ Comparison Report</h1>
        
        <div class="section">
            <h2>1. Read Length Distributions</h2>
            <p>Comparison of read length distributions between real and simulated CCS reads.</p>
            <img src="length_distributions.png" alt="Length Distributions">
        </div>
        
        <div class="section">
            <h2>2. Quality Score Analysis</h2>
            <p>Distribution of mean quality scores per read and their relationship with read length.</p>
            <img src="quality_distributions.png" alt="Quality Distributions">
        </div>
        
        <div class="section">
            <h2>3. GC Content</h2>
            <p>GC content distribution across all reads.</p>
            <img src="gc_distribution.png" alt="GC Distribution">
        </div>
        
        <div class="section">
            <h2>4. Homopolymer Analysis</h2>
            <p>Distribution of homopolymer run lengths (2-15bp) in both datasets.</p>
            <img src="homopolymer_distribution.png" alt="Homopolymer Distribution">
        </div>
        
        <div class="section">
            <h2>5. Summary Statistics</h2>
            <p>See <a href="basic_statistics.tsv">basic_statistics.tsv</a> for detailed metrics.</p>
        </div>
    </body>
    </html>
    """
    
    with open(output_dir / 'comparison_report.html', 'w') as f:
        f.write(html_content)

def main():
    parser = argparse.ArgumentParser(description='Compare real and simulated PacBio CCS FASTQ files')
    parser.add_argument('real_fastq', help='Real PacBio CCS FASTQ file')
    parser.add_argument('sim_fastq', help='Simulated PacBio CCS FASTQ file')
    parser.add_argument('-o', '--output', default='pacbio_comparison', help='Output directory')
    parser.add_argument('-n', '--max-reads', type=int, help='Maximum number of reads to analyze')
    
    args = parser.parse_args()
    
    # Check input files
    if not os.path.exists(args.real_fastq):
        print(f"Error: Real FASTQ file not found: {args.real_fastq}")
        sys.exit(1)
    if not os.path.exists(args.sim_fastq):
        print(f"Error: Simulated FASTQ file not found: {args.sim_fastq}")
        sys.exit(1)
    
    print("Parsing real FASTQ file...")
    real_data = parse_fastq(args.real_fastq, args.max_reads)
    print(f"  Found {len(real_data['lengths'])} reads")
    
    print("Parsing simulated FASTQ file...")
    sim_data = parse_fastq(args.sim_fastq, args.max_reads)
    print(f"  Found {len(sim_data['lengths'])} reads")
    
    print("Creating comparison plots...")
    plot_comparison(real_data, sim_data, args.output)
    
    print("Saving statistics...")
    save_statistics(real_data, sim_data, args.output)
    
    print("Generating report...")
    generate_report(args.output)
    
    print(f"\nAnalysis complete! Results saved to: {args.output}/")
    print(f"Open {args.output}/comparison_report.html to view the report.")

if __name__ == "__main__":
    main()
