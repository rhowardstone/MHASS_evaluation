#!/usr/bin/env python3
"""
Simplified analysis of PacBio HiFi number-of-passes (NP) distribution.

This script performs marginal distribution fitting (NB, Gamma, Log-Normal)
for the number of passes in PacBio HiFi data, with separate panels for different amplicon types.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats, optimize
import argparse
import sys
import warnings

def load_data(counts_file, full_file=None):
    """Load NP count data and optionally full NP/length data."""
    # Load count data
    try:
        df_counts = pd.read_csv(counts_file, sep='\t', header=None, names=['np', 'count'])
        
        # Validate data
        if len(df_counts) < 2:
            raise ValueError("Data file must contain at least 2 rows")
        
        k_values = df_counts['np'].values
        counts = df_counts['count'].values
        
        # Load full data if provided (for statistics only)
        df_full = None
        if full_file:
            df_full = pd.read_csv(full_file, sep='\t')
            print(f"Loaded {len(df_full)} reads with NP and read length data")
            
        return k_values, counts, df_full
        
    except Exception as e:
        print(f"Error loading data: {e}", file=sys.stderr)
        sys.exit(1)

def compute_moments(k_values, counts):
    """Compute mean and variance from count data."""
    total = np.sum(counts)
    mean = np.sum(k_values * counts) / total
    variance = np.sum(counts * (k_values - mean)**2) / total
    return mean, variance

def fit_negative_binomial_mle(k_values, counts):
    """Fit negative binomial using closed-form MLE."""
    k_min = k_values.min()
    k_shifted = k_values - k_min
    
    mean_shifted, variance_shifted = compute_moments(k_shifted, counts)
    
    if variance_shifted <= mean_shifted:
        print("Warning: Variance <= Mean, negative binomial may not be appropriate")
        p = 0.99
    else:
        p = mean_shifted / variance_shifted
    
    r = mean_shifted * p / (1 - p)
    return r, p, k_min

def discretize_continuous(cdf_func, k_values, *params, **kwargs):
    """Properly discretize a continuous distribution."""
    k_min = k_values[0]
    probs = np.zeros(len(k_values))
    
    for i, k in enumerate(k_values):
        if k == k_min:
            probs[i] = cdf_func(k + 0.5, *params, **kwargs)
        else:
            probs[i] = cdf_func(k + 0.5, *params, **kwargs) - cdf_func(k - 0.5, *params, **kwargs)
    
    probs = probs / np.sum(probs)
    return probs

def fit_gamma_mle(k_values, counts):
    """Fit gamma distribution using MLE."""
    mean, variance = compute_moments(k_values, counts)
    scale = variance / mean
    alpha = mean / scale
    
    def neg_log_likelihood(params):
        alpha, scale = params
        if alpha <= 0 or scale <= 0:
            return 1e10
        probs = discretize_continuous(stats.gamma.cdf, k_values, alpha, scale=scale)
        ll = np.sum(counts * np.log(probs + 1e-10))
        return -ll
    
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        result = optimize.minimize(neg_log_likelihood, [alpha, scale], 
                                 bounds=[(0.1, 100), (0.1, 100)],
                                 method='L-BFGS-B')
    return result.x

def fit_lognormal_mle(k_values, counts):
    """Fit log-normal distribution using MLE."""
    total = np.sum(counts)
    log_k = np.log(k_values)
    mu = np.sum(counts * log_k) / total
    sigma = np.sqrt(np.sum(counts * (log_k - mu)**2) / total)
    
    def neg_log_likelihood(params):
        mu, sigma = params
        if sigma <= 0:
            return 1e10
        probs = discretize_continuous(stats.lognorm.cdf, k_values, sigma, scale=np.exp(mu))
        ll = np.sum(counts * np.log(probs + 1e-10))
        return -ll
    
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=RuntimeWarning)
        result = optimize.minimize(neg_log_likelihood, [mu, sigma], 
                                 bounds=[(-10, 10), (0.01, 5)],
                                 method='L-BFGS-B')
    return result.x

def calculate_ks_statistic(k_values, counts, probs):
    """Calculate Kolmogorov-Smirnov statistic."""
    n_total = np.sum(counts)
    # Create empirical CDF
    emp_cdf = np.cumsum(counts) / n_total
    # Create theoretical CDF
    theo_cdf = np.cumsum(probs)
    # KS statistic
    ks_stat = np.max(np.abs(emp_cdf - theo_cdf))
    return ks_stat

def fit_all_models(k_values, counts):
    """Fit all distribution models to the data."""
    fitted_models = {}
    
    # Negative Binomial
    r, p, k_min = fit_negative_binomial_mle(k_values, counts)
    k_shifted = k_values - k_min
    nb_probs = stats.nbinom.pmf(k_shifted, r, p)
    nb_probs = nb_probs / np.sum(nb_probs)
    fitted_models['Negative Binomial'] = ((r, p), nb_probs)
    
    # Gamma
    alpha, scale = fit_gamma_mle(k_values, counts)
    gamma_probs = discretize_continuous(stats.gamma.cdf, k_values, alpha, scale=scale)
    fitted_models['Gamma'] = ((alpha, scale), gamma_probs)
    
    # Log-Normal
    mu, sigma = fit_lognormal_mle(k_values, counts)
    lognorm_probs = discretize_continuous(stats.lognorm.cdf, k_values, sigma, scale=np.exp(mu))
    fitted_models['Log-Normal'] = ((mu, sigma), lognorm_probs)
    
    return fitted_models

def plot_distribution_panel(ax, k_values, counts, fitted_models, title, xlim=60):
    """Plot distribution fits in a single panel."""
    n_total = np.sum(counts)
    
    # Plot observed data
    ax.bar(k_values, counts, width=0.8, color='lightgray', 
           edgecolor='black', linewidth=0.5, label='Observed data')
    
    colors = {'Negative Binomial': 'red', 'Gamma': 'blue', 'Log-Normal': 'green'}
    styles = {'Negative Binomial': '-', 'Gamma': '--', 'Log-Normal': '-.'}
    
    # Calculate KS statistics for all models
    ks_stats = {}
    for name, (params, probs) in fitted_models.items():
        ks_stats[name] = calculate_ks_statistic(k_values, counts, probs)
    
    # Find the best model (lowest KS statistic)
    best_model = min(ks_stats, key=ks_stats.get)
    
    # Plot fitted distributions
    for name, (params, probs) in fitted_models.items():
        expected_counts = probs * n_total
        ks_stat = ks_stats[name]
        
        # Create label with parameters and KS statistic
        if name == 'Negative Binomial':
            r, p = params
            label = f'{name} (r={r:.2f}, p={p:.3f})'
        elif name == 'Gamma':
            alpha, scale = params
            label = f'{name} (α={alpha:.2f}, θ={scale:.2f})'
        elif name == 'Log-Normal':
            mu, sigma = params
            label = f'{name} (μ={mu:.2f}, σ={sigma:.2f})'
        
        # Add KS statistic to label, with star for best model
        if name == best_model:
            label += f' | KS={ks_stat:.4f} ★'
        else:
            label += f' | KS={ks_stat:.4f}'
        
        # Plot with thicker line for best model
        linewidth = 3.5 if name == best_model else 2.5
        ax.plot(k_values, expected_counts, color=colors[name], 
                linestyle=styles[name], linewidth=linewidth, label=label)
    
    ax.set_xlabel('Number of Passes (NP)', fontsize=14, fontweight='bold')
    ax.set_ylabel('Count', fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold')
    
    # Customize legend
    legend = ax.legend(loc='upper right', fontsize=10, frameon=True, fancybox=True, shadow=True)
    
    # Make legend text bold
    for text in legend.get_texts():
        text.set_weight('bold')
    
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, xlim)
    
    # Make tick labels bold
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_weight('bold')

def print_summary_for_dataset(name, k_values, counts, fitted_models):
    """Print summary for a single dataset."""
    n_total = np.sum(counts)
    mean, variance = compute_moments(k_values, counts)
    
    print(f"\n{name} Data Summary:")
    print(f"  Total reads: {n_total:,}")
    print(f"  Mean NP: {mean:.2f}")
    print(f"  Variance: {variance:.2f}")
    print(f"  Std Dev: {np.sqrt(variance):.2f}")
    print(f"  Range: {k_values.min()} - {k_values.max()}")
    
    # Calculate KS stats and find best model
    model_stats = []
    for model_name, (params, probs) in fitted_models.items():
        ks_stat = calculate_ks_statistic(k_values, counts, probs)
        if model_name == 'Log-Normal':
            model_stats.append((model_name, params, ks_stat))
    
    # Print log-normal parameters
    for model_name, params, ks_stat in model_stats:
        if model_name == 'Log-Normal':
            print(f"\n{name} Log-Normal parameters for MHASS:")
            print(f"  --lognormal_mu {params[0]:.2f}")
            print(f"  --lognormal_sigma {params[1]:.2f}")

def main():
    parser = argparse.ArgumentParser(
        description='Analysis of PacBio HiFi NP distribution for multiple amplicon types'
    )
    parser.add_argument('--titan_counts', required=True, help='TSV file with Titan NP counts')
    parser.add_argument('--s16_counts', required=True, help='TSV file with 16S NP counts')
    parser.add_argument('--titan_full', help='Full data file for Titan (optional)')
    parser.add_argument('--s16_full', help='Full data file for 16S (optional)')
    parser.add_argument('--out', default='np_distribution_comparison.png',
                       help='Output plot filename')
    
    args = parser.parse_args()
    
    # Set up bold fonts globally
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titleweight'] = 'bold'
    plt.rcParams['legend.fontsize'] = 10
    
    # Load data for both amplicon types
    print("Loading Titan amplicon data...")
    titan_k, titan_counts, titan_full = load_data(args.titan_counts, args.titan_full)
    
    print("Loading 16S amplicon data...")
    s16_k, s16_counts, s16_full = load_data(args.s16_counts, args.s16_full)
    
    # Fit models for both datasets
    print("\nFitting distributions for Titan amplicons...")
    titan_models = fit_all_models(titan_k, titan_counts)
    
    print("Fitting distributions for 16S amplicons...")
    s16_models = fit_all_models(s16_k, s16_counts)
    
    # Create 2x1 plot
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 14))
    
    # Plot Titan data (top panel)
    plot_distribution_panel(ax1, titan_k, titan_counts, titan_models, 
                           'Titan Amplicons (2000-3000 bp)', xlim=60)
    
    # Plot 16S data (bottom panel)
    plot_distribution_panel(ax2, s16_k, s16_counts, s16_models, 
                           '16S Amplicons (1400-1700 bp)', xlim=60)
    
    plt.tight_layout()
    plt.savefig(args.out, dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print summaries
    print("\n" + "="*80)
    print("DISTRIBUTION ANALYSIS SUMMARY")
    print("="*80)
    
    print_summary_for_dataset('Titan', titan_k, titan_counts, titan_models)
    print_summary_for_dataset('16S', s16_k, s16_counts, s16_models)
    
    print(f"\nResults saved to {args.out}")
    print("Done!")

if __name__ == "__main__":
    main()
