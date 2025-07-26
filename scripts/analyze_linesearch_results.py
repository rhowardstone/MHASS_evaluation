#!/usr/bin/env python3
"""
analyze_linesearch_results.py
Analyze line search results for subread accuracy optimization.
Parallelized version to process multiple parameter combinations simultaneously.
"""

import argparse
import os
import glob
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import subprocess
import sys
from multiprocessing import Pool, cpu_count
import time

def parse_args():
    parser = argparse.ArgumentParser(
        description='Analyze line search results for subread accuracy optimization'
    )
    parser.add_argument('--base_dir', required=True,
                        help='Base directory containing all Analysis-* folders')
    parser.add_argument('--real_closest', required=True,
                        help='Path to real data closest file')
    parser.add_argument('--output_dir', required=True,
                        help='Output directory for analysis results')
    parser.add_argument('--measure_script', required=True,
                        help='Path to measure_distribution_differences.py')
    parser.add_argument('--skip_measurement', action='store_true',
                        help='Skip measurement step if already done')
    parser.add_argument('--threads', type=int, default=None,
                        help='Number of parallel threads (default: auto-detect)')
    return parser.parse_args()

def run_distribution_measurement(measure_script, real_closest, sim_closest, output_prefix):
    """Run the distribution measurement script."""
    cmd = [
        sys.executable, measure_script,
        '--real_closest', real_closest,
        '--sim_closest', sim_closest,
        '--output_prefix', output_prefix,
        '--plot'
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return (True, output_prefix, None)
    except subprocess.CalledProcessError as e:
        return (False, output_prefix, f"Error: {e}\nstderr: {e.stderr}")

def extract_params_from_dirname(dirname):
    """Extract parameters from directory name like Analysis-acc0.75_mu3.08_sigma0.94"""
    parts = dirname.split('-')[-1].split('_')
    params = {}
    
    for part in parts:
        if part.startswith('acc'):
            params['accuracy'] = float(part[3:])
        elif part.startswith('mu'):
            params['mu'] = float(part[2:])
        elif part.startswith('sigma'):
            params['sigma'] = float(part[5:])
    
    return params

def process_single_directory(args):
    """Process a single analysis directory. Used for parallel processing."""
    analysis_dir, real_closest, output_dir, measure_script, skip_measurement = args
    
    dirname = os.path.basename(analysis_dir)
    
    # Extract parameters
    params = extract_params_from_dirname(dirname)
    if not params:
        return None
    
    # Check if simulation was successful
    sim_closest = os.path.join(analysis_dir, "simulated/ref-closest-sim.txt")
    if not os.path.exists(sim_closest):
        return None
    
    # Measure distribution differences
    metrics_file = os.path.join(output_dir, f"{dirname}_metrics.json")
    
    if not skip_measurement or not os.path.exists(metrics_file):
        output_prefix = os.path.join(output_dir, dirname)
        success, _, error_msg = run_distribution_measurement(
            measure_script, real_closest, sim_closest, output_prefix
        )
        
        if not success:
            print(f"Failed to measure {dirname}: {error_msg}")
            return None
    
    # Load metrics
    if os.path.exists(metrics_file):
        with open(metrics_file, 'r') as f:
            metrics = json.load(f)
        
        # Combine parameters and metrics
        result = {**params, **metrics}
        result['directory'] = analysis_dir
        return result
    else:
        return None

def collect_all_results_parallel(base_dir, real_closest, output_dir, measure_script, 
                                skip_measurement=False, n_threads=None):
    """Collect results from all line search runs using parallel processing."""
    
    # Find all Analysis-* directories with the new naming convention
    analysis_dirs = sorted(glob.glob(os.path.join(base_dir, "Analysis-acc*_mu*_sigma*")))
    
    print(f"Found {len(analysis_dirs)} analysis directories")
    
    # Determine number of threads
    if n_threads is None:
        n_threads = min(cpu_count(), len(analysis_dirs))
    else:
        n_threads = min(n_threads, len(analysis_dirs))
    
    print(f"Using {n_threads} parallel threads")
    
    # Prepare arguments for parallel processing
    process_args = [
        (analysis_dir, real_closest, output_dir, measure_script, skip_measurement)
        for analysis_dir in analysis_dirs
    ]
    
    # Process in parallel
    start_time = time.time()
    results = []
    
    with Pool(n_threads) as pool:
        # Use imap_unordered for better progress tracking
        total = len(process_args)
        completed = 0
        
        for result in pool.imap_unordered(process_single_directory, process_args):
            completed += 1
            if completed % 5 == 0 or completed == total:
                elapsed = time.time() - start_time
                rate = completed / elapsed
                remaining = (total - completed) / rate if rate > 0 else 0
                print(f"Progress: {completed}/{total} ({completed/total*100:.1f}%) - "
                      f"Est. remaining: {remaining/60:.1f} min")
            
            if result is not None:
                results.append(result)
    
    elapsed_total = time.time() - start_time
    print(f"\nCompleted processing in {elapsed_total/60:.1f} minutes")
    print(f"Successfully processed {len(results)} out of {len(analysis_dirs)} directories")
    
    return pd.DataFrame(results)

def create_line_plots(df, output_dir):
    """Create line plots showing how metrics change with subread accuracy."""
    metrics_to_plot = [
        ('kl_symmetric', 'Symmetric KL Divergence'),
        ('js_divergence', 'Jensen-Shannon Divergence'),
        ('wasserstein_distance', 'Wasserstein Distance'),
        ('total_variation', 'Total Variation Distance'),
        ('mean_diff', 'Mean Difference'),
        ('std_diff', 'Standard Deviation Difference')
    ]
    
    # Create a figure with multiple subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    # Sort by accuracy for proper line plotting
    df_sorted = df.sort_values('accuracy')
    
    for i, (metric, title) in enumerate(metrics_to_plot):
        if metric not in df.columns:
            continue
        
        ax = axes[i]
        
        # Plot the main line
        ax.plot(df_sorted['accuracy'], df_sorted[metric], 'b-o', linewidth=2, 
                markersize=6, label=metric)
        
        # Highlight the minimum
        min_idx = df_sorted[metric].idxmin()
        min_row = df_sorted.loc[min_idx]
        ax.plot(min_row['accuracy'], min_row[metric], 'ro', markersize=10, 
                label=f'Best: {min_row["accuracy"]:.2f}')
        
        # Add text annotation for best value
        ax.annotate(f'Best: {min_row["accuracy"]:.2f}\nValue: {min_row[metric]:.4f}',
                   xy=(min_row['accuracy'], min_row[metric]),
                   xytext=(10, 10), textcoords='offset points',
                   bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.7),
                   arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        
        # Customize plot
        ax.set_xlabel('Subread Accuracy')
        ax.set_ylabel(title)
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend()
        
        # Set x-axis to show all accuracy values
        ax.set_xticks(np.arange(0.05, 1.05, 0.1))
        ax.set_xlim(0.0, 1.05)
    
    plt.suptitle('Effect of Subread Accuracy on Distribution Metrics', fontsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'subread_accuracy_effects.png'), 
               dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create a combined plot showing all metrics normalized
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Normalize each metric to 0-1 range for comparison
    colors = plt.cm.Set1(np.linspace(0, 1, len(metrics_to_plot)))
    
    for i, (metric, title) in enumerate(metrics_to_plot):
        if metric not in df.columns:
            continue
        
        # Normalize metric values
        metric_values = df_sorted[metric].values
        normalized_values = (metric_values - metric_values.min()) / (metric_values.max() - metric_values.min())
        
        ax.plot(df_sorted['accuracy'], normalized_values, 'o-', 
               color=colors[i], linewidth=2, markersize=6, label=title)
    
    ax.set_xlabel('Subread Accuracy')
    ax.set_ylabel('Normalized Metric Value (0=best, 1=worst)')
    ax.set_title('Normalized Comparison of All Metrics vs Subread Accuracy')
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.set_xticks(np.arange(0.05, 1.05, 0.1))
    ax.set_xlim(0.0, 1.05)
    ax.set_ylim(-0.05, 1.05)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'normalized_metrics_comparison.png'), 
               dpi=300, bbox_inches='tight')
    plt.close()

def find_best_parameters(df, metric='kl_symmetric'):
    """Find the best parameter combination based on a metric."""
    # Lower is better for distance metrics
    best_idx = df[metric].idxmin()
    best_params = df.loc[best_idx]
    
    return best_params

def create_summary_report(df, output_dir):
    """Create a summary report of the line search results."""
    report_file = os.path.join(output_dir, 'linesearch_summary_report.txt')
    
    with open(report_file, 'w') as f:
        f.write("Line Search Summary Report\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Total subread accuracy values tested: {len(df)}\n")
        f.write(f"Accuracy range: {df['accuracy'].min():.2f} - {df['accuracy'].max():.2f}\n")
        f.write(f"Fixed lognormal parameters: mu={df['mu'].iloc[0]:.2f}, sigma={df['sigma'].iloc[0]:.2f}\n\n")
        
        # Best parameters for each metric
        f.write("Best Subread Accuracy by Metric:\n")
        f.write("-" * 40 + "\n")
        
        metrics = ['kl_symmetric', 'js_divergence', 'wasserstein_distance', 'total_variation']
        for metric in metrics:
            if metric not in df.columns:
                continue
                
            best = find_best_parameters(df, metric)
            f.write(f"\n{metric.replace('_', ' ').title()}:\n")
            f.write(f"  Best accuracy: {best['accuracy']:.3f}\n")
            f.write(f"  Metric value: {best[metric]:.6f}\n")
        
        # Statistics across all accuracy values
        f.write("\n\nMetric Statistics Across All Accuracy Values:\n")
        f.write("-" * 50 + "\n")
        
        for metric in metrics:
            if metric not in df.columns:
                continue
            f.write(f"\n{metric.replace('_', ' ').title()}:\n")
            f.write(f"  Mean: {df[metric].mean():.6f}\n")
            f.write(f"  Std:  {df[metric].std():.6f}\n")
            f.write(f"  Min:  {df[metric].min():.6f} (at accuracy {df.loc[df[metric].idxmin(), 'accuracy']:.3f})\n")
            f.write(f"  Max:  {df[metric].max():.6f} (at accuracy {df.loc[df[metric].idxmax(), 'accuracy']:.3f})\n")
        
        # Top 5 best combinations
        f.write("\n\nTop 5 Subread Accuracy Values (by symmetric KL divergence):\n")
        f.write("-" * 60 + "\n")
        top5 = df.nsmallest(5, 'kl_symmetric')[['accuracy', 'kl_symmetric', 
                                                'js_divergence', 'wasserstein_distance']]
        f.write(top5.to_string(index=False))
        
        # Worst 5 combinations for comparison
        f.write("\n\nWorst 5 Subread Accuracy Values (by symmetric KL divergence):\n")
        f.write("-" * 60 + "\n")
        worst5 = df.nlargest(5, 'kl_symmetric')[['accuracy', 'kl_symmetric', 
                                                 'js_divergence', 'wasserstein_distance']]
        f.write(worst5.to_string(index=False))
        
    print(f"\nSummary report saved to: {report_file}")

def main():
    args = parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    print("Collecting and analyzing line search results...")
    
    # Collect all results with parallel processing
    df = collect_all_results_parallel(
        args.base_dir, args.real_closest, args.output_dir, 
        args.measure_script, args.skip_measurement, args.threads
    )
    
    if len(df) == 0:
        print("No results found!")
        return
    
    print(f"\nCollected results from {len(df)} parameter combinations")
    
    # Save raw results
    results_file = os.path.join(args.output_dir, 'all_linesearch_results.csv')
    df.to_csv(results_file, index=False)
    print(f"Raw results saved to: {results_file}")
    
    # Create visualizations
    print("\nCreating line plots...")
    create_line_plots(df, args.output_dir)
    
    # Find best parameters
    print("\nFinding best parameters...")
    best_params = find_best_parameters(df)
    print(f"Best subread accuracy (by symmetric KL divergence): {best_params['accuracy']:.3f}")
    print(f"KL divergence: {best_params['kl_symmetric']:.6f}")
    
    # Create summary report
    create_summary_report(df, args.output_dir)
    
    print(f"\nAnalysis complete! Results saved to: {args.output_dir}")

if __name__ == '__main__':
    main()
