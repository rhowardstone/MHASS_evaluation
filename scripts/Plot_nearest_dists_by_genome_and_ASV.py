#!/usr/bin/env python3
"""
Plot_nearest_dists_by_genome_and_ASV.py
Create per-genome KDE plots of edit-distance-to-nearest-ground-truth ASV.

Outputs
  <genome>_real.png   – real reads vs. consensus ASVs
  <genome>_sim.png    – simulated reads vs. reference ASVs

Legend entry shows the actual ASV ID (n=<multiplicity>)
  e.g. Escherichia_coli-1-2-3 (n=87)
"""

import argparse, os, re, gzip
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def extract_genome(seq_id):
    """Extract genome from deduped ID format like 'Escherichia_coli-1-2-3'."""
    if pd.isna(seq_id) or not isinstance(seq_id, str):
        return None
    
    # Find the position of the first dash followed by a digit
    match = re.search(r'^(.+?)-\d', seq_id)
    if match:
        return match.group(1)
    return seq_id

def extract_indices(seq_id):
    """Extract indices from deduped ID format like 'Escherichia_coli-1-2-3'."""
    if pd.isna(seq_id) or not isinstance(seq_id, str):
        return None
    
    # Extract all the numeric indices after the genome name
    match = re.search(r'-(\d+(?:-\d+)*)$', seq_id)
    if match:
        return match.group(1)
    return None

def load_closest(fname: str) -> pd.DataFrame:
    rows = []
    opener = gzip.open if fname.endswith((".gz", ".gzip")) else open
    with opener(fname, 'rt') as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            p = ln.split('\t')
            if len(p) < 3:
                continue
            read, dist = p[0], int(p[1])
            for ref in p[2:]:
                rows.append((read, ref, dist))
    return pd.DataFrame(rows, columns=["read_id", "ref_id", "dist"])

def prep_data(df):
    """Prepare data with clean ID format."""
    df = df.copy()
    
    # Extract genome from clean ID
    df["genome"] = df["ref_id"].apply(extract_genome)
    
    # Count multiplicity by reference
    df["multiplicity"] = df.groupby("ref_id")["ref_id"].transform("count")
    
    # Drop rows with missing values
    return df.dropna(subset=["genome"])

def plot_one(genome: str, data: pd.DataFrame, kind: str, outdir: str):
    sns.set_style("whitegrid")
    plt.figure(figsize=(10, 6))

    total_reads = data["multiplicity"].sum()

    # Sort ASVs by their ID for consistent ordering
    sorted_asvs = sorted(data["ref_id"].unique())
    
    for ref_id in sorted_asvs:
        sub = data[data["ref_id"] == ref_id]
        mult = sub["multiplicity"].iloc[0]
        weight = mult / total_reads  # mixture weight
        
        # Use the actual ASV ID in the label
        label = f"{ref_id} (n={mult})"

        # Use KDE plot with proper parameters
        distances = sub["dist"].to_numpy()
        weights_array = [weight] * len(sub)
        
        if len(distances) > 1:  # Need at least 2 points for KDE
            sns.kdeplot(
                x=distances,
                weights=weights_array,
                bw_adjust=0.5,
                common_norm=False,
                label=label,
            )
        else:
            # For single points, just plot a vertical line
            plt.axvline(x=distances[0], label=label, alpha=0.7)

    plt.title(f"{genome} – {kind}")
    plt.xlabel("Edit distance to nearest ground-truth ASV")
    plt.ylabel("Density")
    plt.legend(fontsize="small", frameon=False)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{genome}_{kind}.png"), dpi=300)
    plt.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real_closest", required=True)
    ap.add_argument("--sim_closest",  required=True)
    ap.add_argument("--outdir",       required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading closest distance files...")
    df_real = load_closest(args.real_closest)
    df_sim  = load_closest(args.sim_closest)
    
    print(f"Real data: {len(df_real)} entries")
    print(f"Simulated data: {len(df_sim)} entries")
    
    print("Preparing data...")
    df_real = prep_data(df_real)
    df_sim  = prep_data(df_sim)
    
    print(f"After preparation - Real: {len(df_real)} entries, Sim: {len(df_sim)} entries")
    
    # Get all unique genomes
    all_genomes = sorted(set(df_real["genome"].unique()) | set(df_sim["genome"].unique()))
    print(f"Found {len(all_genomes)} genomes: {', '.join(all_genomes)}")
    
    for g in all_genomes:
        real_subset = df_real[df_real["genome"] == g]
        sim_subset = df_sim[df_sim["genome"] == g]
        
        if len(real_subset) > 0:
            print(f"Plotting {g}_real.png ({len(real_subset)} reads)")
            plot_one(g, real_subset, "real", args.outdir)
        
        if len(sim_subset) > 0:
            print(f"Plotting {g}_sim.png ({len(sim_subset)} reads)")
            plot_one(g, sim_subset, "sim", args.outdir)

    print("✅ Plots written to", args.outdir)


if __name__ == "__main__":
    main()
