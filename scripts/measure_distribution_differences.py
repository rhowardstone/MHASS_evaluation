#!/usr/bin/env python3
"""
measure_distribution_differences.py
Measure distribution differences between simulated and real data edit distances.
Calculates multiple metrics including KL divergence, JS divergence, Wasserstein distance, etc.
"""

import argparse
import os
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import jensenshannon
from scipy.stats import wasserstein_distance, energy_distance
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import json


def parse_args():
    parser = argparse.ArgumentParser(
        description='Measure distribution differences between real and simulated edit distances'
    )
    parser.add_argument('--real_closest', required=True,
                        help='Path to real data closest file (consensus-closest-real.txt)')
    parser.add_argument('--sim_closest', required=True,
                        help='Path to simulated data closest file (ref-closest-sim.txt)')
    parser.add_argument('--output_prefix', required=True,
                        help='Prefix for output files')
    parser.add_argument('--max_distance', type=int, default=None,
                        help='Maximum edit distance to consider')
    parser.add_argument('--bins', type=int, default=100,
                        help='Number of bins for histogram comparison (default: 100)')
    parser.add_argument('--plot', action='store_true',
                        help='Generate comparison plots')
    return parser.parse_args()


def read_closest_file(filepath, max_distance=None):
    """Read closest file and extract edit distances."""
    distances = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                distance = int(parts[1])
                if max_distance is None or distance <= max_distance:
                    distances.append(distance)
    return np.array(distances)


def create_probability_distributions(real_dists, sim_dists, bins=100):
    """Create normalized histograms (probability distributions) for comparison."""
    # Determine the range for binning
    min_val = min(real_dists.min(), sim_dists.min())
    max_val = max(real_dists.max(), sim_dists.max())
    
    # Create bins
    bin_edges = np.linspace(min_val, max_val + 1, bins + 1)
    
    # Create histograms
    real_hist, _ = np.histogram(real_dists, bins=bin_edges)
    sim_hist, _ = np.histogram(sim_dists, bins=bin_edges)
    
    # Normalize to create probability distributions
    real_prob = real_hist / real_hist.sum()
    sim_prob = sim_hist / sim_hist.sum()
    
    return real_prob, sim_prob, bin_edges


def kl_divergence(p, q, epsilon=1e-10):
    """Calculate KL divergence KL(P||Q)."""
    # Add epsilon to avoid log(0)
    p = p + epsilon
    q = q + epsilon
    
    # Renormalize
    p = p / p.sum()
    q = q / q.sum()
    
    # Calculate KL divergence
    kl = np.sum(p * np.log(p / q))
    return kl


def calculate_distribution_metrics(real_dists, sim_dists, bins=100):
    """Calculate various distribution comparison metrics."""
    metrics = {}
    
    # Basic statistics
    metrics['real_mean'] = float(np.mean(real_dists))
    metrics['real_std'] = float(np.std(real_dists))
    metrics['real_median'] = float(np.median(real_dists))
    metrics['real_q25'] = float(np.percentile(real_dists, 25))
    metrics['real_q75'] = float(np.percentile(real_dists, 75))
    
    metrics['sim_mean'] = float(np.mean(sim_dists))
    metrics['sim_std'] = float(np.std(sim_dists))
    metrics['sim_median'] = float(np.median(sim_dists))
    metrics['sim_q25'] = float(np.percentile(sim_dists, 25))
    metrics['sim_q75'] = float(np.percentile(sim_dists, 75))
    
    # Differences in basic stats
    metrics['mean_diff'] = abs(metrics['real_mean'] - metrics['sim_mean'])
    metrics['std_diff'] = abs(metrics['real_std'] - metrics['sim_std'])
    metrics['median_diff'] = abs(metrics['real_median'] - metrics['sim_median'])
    
    # Create probability distributions
    real_prob, sim_prob, bin_edges = create_probability_distributions(real_dists, sim_dists, bins)
    
    # KL Divergence (both directions)
    metrics['kl_real_to_sim'] = float(kl_divergence(real_prob, sim_prob))
    metrics['kl_sim_to_real'] = float(kl_divergence(sim_prob, real_prob))
    metrics['kl_symmetric'] = (metrics['kl_real_to_sim'] + metrics['kl_sim_to_real']) / 2
    
    # Jensen-Shannon Divergence (symmetric version of KL)
    metrics['js_divergence'] = float(jensenshannon(real_prob, sim_prob) ** 2)  # Squared JS distance
    
    # Wasserstein Distance (Earth Mover's Distance)
    metrics['wasserstein_distance'] = float(wasserstein_distance(real_dists, sim_dists))
    
    # Energy Distance
    metrics['energy_distance'] = float(energy_distance(real_dists, sim_dists))
    
    # Kolmogorov-Smirnov test
    ks_stat, ks_pvalue = stats.ks_2samp(real_dists, sim_dists)
    metrics['ks_statistic'] = float(ks_stat)
    metrics['ks_pvalue'] = float(ks_pvalue)
    
    # Skip chi-square test - it's not well-suited for comparing distributions 
    # with different sample sizes. We have many other better metrics.
    metrics['chi2_statistic'] = np.nan
    metrics['chi2_pvalue'] = np.nan
    
    # Total Variation Distance
    metrics['total_variation'] = float(0.5 * np.sum(np.abs(real_prob - sim_prob)))
    
    # Hellinger Distance
    metrics['hellinger_distance'] = float(np.sqrt(0.5 * np.sum((np.sqrt(real_prob) - np.sqrt(sim_prob))**2)))
    
    # Bhattacharyya Distance
    bc_coefficient = np.sum(np.sqrt(real_prob * sim_prob))
    metrics['bhattacharyya_distance'] = float(-np.log(bc_coefficient + 1e-10))
    
    return metrics


def create_comparison_plots(real_dists, sim_dists, output_prefix):
    """Create comparison plots between real and simulated distributions."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Overlaid histograms
    ax = axes[0, 0]
    bins = np.linspace(0, max(real_dists.max(), sim_dists.max()), 50)
    ax.hist(real_dists, bins=bins, alpha=0.5, label='Real', density=True, color='blue')
    ax.hist(sim_dists, bins=bins, alpha=0.5, label='Simulated', density=True, color='red')
    ax.set_xlabel('Edit Distance')
    ax.set_ylabel('Density')
    ax.set_title('Distribution Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. CDFs
    ax = axes[0, 1]
    sorted_real = np.sort(real_dists)
    sorted_sim = np.sort(sim_dists)
    ax.plot(sorted_real, np.arange(len(sorted_real)) / len(sorted_real), 
            label='Real', color='blue', linewidth=2)
    ax.plot(sorted_sim, np.arange(len(sorted_sim)) / len(sorted_sim), 
            label='Simulated', color='red', linewidth=2)
    ax.set_xlabel('Edit Distance')
    ax.set_ylabel('Cumulative Probability')
    ax.set_title('Cumulative Distribution Functions')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 3. Q-Q plot
    ax = axes[1, 0]
    stats.probplot(sim_dists, dist=stats.norm, sparams=(real_dists.mean(), real_dists.std()), 
                   plot=ax)
    ax.set_title('Q-Q Plot (Simulated vs Real)')
    ax.grid(True, alpha=0.3)
    
    # 4. Box plots
    ax = axes[1, 1]
    data_to_plot = [real_dists, sim_dists]
    box_plot = ax.boxplot(data_to_plot, labels=['Real', 'Simulated'], patch_artist=True)
    box_plot['boxes'][0].set_facecolor('blue')
    box_plot['boxes'][1].set_facecolor('red')
    ax.set_ylabel('Edit Distance')
    ax.set_title('Distribution Summary')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_comparison_plots.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Additional plot: KDE comparison
    plt.figure(figsize=(10, 6))
    sns.kdeplot(real_dists, label='Real', color='blue', fill=True, alpha=0.3)
    sns.kdeplot(sim_dists, label='Simulated', color='red', fill=True, alpha=0.3)
    plt.xlabel('Edit Distance')
    plt.ylabel('Density')
    plt.title('Kernel Density Estimation Comparison')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{output_prefix}_kde_comparison.png", dpi=300, bbox_inches='tight')
    plt.close()


def main():
    args = parse_args()
    
    # Read distance data
    print(f"Reading real data from: {args.real_closest}")
    real_distances = read_closest_file(args.real_closest, args.max_distance)
    print(f"  Found {len(real_distances)} distances")
    
    print(f"Reading simulated data from: {args.sim_closest}")
    sim_distances = read_closest_file(args.sim_closest, args.max_distance)
    print(f"  Found {len(sim_distances)} distances")
    
    # Calculate metrics
    print("\nCalculating distribution comparison metrics...")
    metrics = calculate_distribution_metrics(real_distances, sim_distances, args.bins)
    
    # Print summary
    print("\n=== Distribution Comparison Summary ===")
    print(f"Real data: mean={metrics['real_mean']:.2f}, std={metrics['real_std']:.2f}, "
          f"median={metrics['real_median']:.2f}")
    print(f"Sim data:  mean={metrics['sim_mean']:.2f}, std={metrics['sim_std']:.2f}, "
          f"median={metrics['sim_median']:.2f}")
    print(f"\nKey metrics:")
    print(f"  KL divergence (real→sim): {metrics['kl_real_to_sim']:.4f}")
    print(f"  KL divergence (sim→real): {metrics['kl_sim_to_real']:.4f}")
    print(f"  JS divergence: {metrics['js_divergence']:.4f}")
    print(f"  Wasserstein distance: {metrics['wasserstein_distance']:.2f}")
    print(f"  Total variation: {metrics['total_variation']:.4f}")
    print(f"  KS statistic: {metrics['ks_statistic']:.4f} (p={metrics['ks_pvalue']:.4e})")
    
    # Save metrics to JSON
    json_file = f"{args.output_prefix}_metrics.json"
    with open(json_file, 'w') as f:
        json.dump(metrics, f, indent=2)
    print(f"\nMetrics saved to: {json_file}")
    
    # Save metrics to CSV for easy reading
    csv_file = f"{args.output_prefix}_metrics.csv"
    pd.DataFrame([metrics]).to_csv(csv_file, index=False)
    print(f"Metrics also saved to: {csv_file}")
    
    # Create plots if requested
    if args.plot:
        print("\nCreating comparison plots...")
        create_comparison_plots(real_distances, sim_distances, args.output_prefix)
        print(f"Plots saved with prefix: {args.output_prefix}")


if __name__ == '__main__':
    main()
