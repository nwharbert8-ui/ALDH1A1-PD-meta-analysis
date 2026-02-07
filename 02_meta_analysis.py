#!/usr/bin/env python3
"""
02_meta_analysis.py
====================
Perform random-effects meta-analysis (DerSimonian-Laird) of Pearson
correlation coefficients using Fisher's z-transformation.

Reads: results/per_dataset_correlations.csv
Produces:
  results/meta_analysis_pooled.csv
  results/meta_analysis_forest_data.csv
  figures/forest_plot_ALDH1A1_TH.pdf (and other key pairs)
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')

# Import shared configuration
import sys
sys.path.insert(0, os.path.dirname(__file__))
from config import get_results_dir, get_figures_dir

RESULTS_DIR = get_results_dir()
FIGURES_DIR = get_figures_dir()


# ---------------------------------------------------------------------------
# Fisher's z-transformation
# ---------------------------------------------------------------------------
def fisher_z(r):
    """Transform Pearson r to Fisher's z."""
    r = np.clip(r, -0.9999, 0.9999)
    return 0.5 * np.log((1 + r) / (1 - r))


def fisher_z_inv(z):
    """Back-transform Fisher's z to Pearson r."""
    return np.tanh(z)


def fisher_z_se(n):
    """Standard error of Fisher's z."""
    return 1.0 / np.sqrt(n - 3)


# ---------------------------------------------------------------------------
# DerSimonian-Laird Random-Effects Meta-Analysis
# ---------------------------------------------------------------------------
def random_effects_meta(z_values, se_values):
    """
    DerSimonian-Laird random-effects meta-analysis.
    
    Parameters
    ----------
    z_values : array-like, Fisher's z-transformed correlations
    se_values : array-like, standard errors of z
    
    Returns
    -------
    dict with: pooled_z, pooled_z_se, pooled_r, ci_lower_r, ci_upper_r,
               Q, I2, tau2, p_heterogeneity, k
    """
    z = np.array(z_values, dtype=float)
    se = np.array(se_values, dtype=float)
    k = len(z)
    
    if k == 0:
        return None
    
    # Fixed-effect weights
    w = 1.0 / (se ** 2)
    
    # Fixed-effect pooled estimate
    z_fe = np.sum(w * z) / np.sum(w)
    
    # Cochran's Q
    Q = np.sum(w * (z - z_fe) ** 2)
    df = k - 1
    p_het = 1 - stats.chi2.cdf(Q, df) if df > 0 else 1.0
    
    # DerSimonian-Laird tau^2
    C = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0, (Q - df) / C) if C > 0 else 0
    
    # I^2
    I2 = max(0, (Q - df) / Q * 100) if Q > 0 else 0
    
    # Random-effects weights
    w_re = 1.0 / (se ** 2 + tau2)
    
    # Pooled estimate
    pooled_z = np.sum(w_re * z) / np.sum(w_re)
    pooled_se = 1.0 / np.sqrt(np.sum(w_re))
    
    # 95% CI
    ci_lower_z = pooled_z - 1.96 * pooled_se
    ci_upper_z = pooled_z + 1.96 * pooled_se
    
    return {
        'pooled_z': pooled_z,
        'pooled_z_se': pooled_se,
        'pooled_r': fisher_z_inv(pooled_z),
        'ci_lower_r': fisher_z_inv(ci_lower_z),
        'ci_upper_r': fisher_z_inv(ci_upper_z),
        'Q': Q,
        'I2': I2,
        'tau2': tau2,
        'p_heterogeneity': p_het,
        'k': k,
    }


# ---------------------------------------------------------------------------
# Forest Plot
# ---------------------------------------------------------------------------
def make_forest_plot(pair_name, forest_data, condition, outpath):
    """Generate a forest plot for a single gene pair in one condition."""
    fig, ax = plt.subplots(figsize=(10, 4 + 0.4 * len(forest_data)))
    
    datasets = forest_data['dataset'].values
    r_vals = forest_data['r'].values
    ci_low = forest_data['ci_lower_r'].values
    ci_high = forest_data['ci_upper_r'].values
    n_vals = forest_data['n'].values
    
    y_pos = np.arange(len(datasets))
    
    # Individual studies
    for i in range(len(datasets)):
        ax.plot([ci_low[i], ci_high[i]], [y_pos[i], y_pos[i]], 'b-', linewidth=1)
        ax.plot(r_vals[i], y_pos[i], 'bs', markersize=6 + n_vals[i] / 5)
    
    # Pooled estimate
    pooled = forest_data['pooled_r'].iloc[0] if 'pooled_r' in forest_data.columns else None
    if pooled is not None:
        ax.axvline(pooled, color='red', linestyle='--', alpha=0.7, label=f'Pooled r = {pooled:.3f}')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{d} (n={n})" for d, n in zip(datasets, n_vals)])
    ax.set_xlabel('Pearson r')
    ax.set_title(f'{pair_name} — {condition}')
    ax.legend(loc='best')
    ax.invert_yaxis()
    
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 80)
    print("ALDH1A1-PD META-ANALYSIS: Script 02 — Random-Effects Meta-Analysis")
    print("=" * 80)
    
    # Load per-dataset correlations
    corr_path = os.path.join(RESULTS_DIR, 'per_dataset_correlations.csv')
    if not os.path.exists(corr_path):
        print("ERROR: per_dataset_correlations.csv not found. Run 01_download_and_correlate.py first.")
        return
    
    df = pd.read_csv(corr_path)
    
    # Get unique target gene pairs
    target_pairs = df[df['category'].isin(['ALDH1A1-DA', 'DA-DA', 'SNCA'])]['pair_name'].unique()
    
    pooled_results = []
    forest_records = []
    
    for pair in sorted(target_pairs):
        pair_data = df[df['pair_name'] == pair].copy()
        category = pair_data['category'].iloc[0]
        
        for condition, r_col, n_col in [('control', 'r_ctrl', 'n_ctrl'),
                                         ('PD', 'r_pd', 'n_pd')]:
            sub = pair_data.dropna(subset=[r_col, n_col])
            sub = sub[sub[n_col] >= 5]  # Minimum sample size
            
            if len(sub) == 0:
                continue
            
            r_values = sub[r_col].values
            n_values = sub[n_col].values.astype(int)
            
            # Fisher's z-transform
            z_values = fisher_z(r_values)
            se_values = fisher_z_se(n_values)
            
            # Meta-analysis
            meta = random_effects_meta(z_values, se_values)
            
            if meta is None:
                continue
            
            pooled_results.append({
                'pair_name': pair,
                'category': category,
                'condition': condition,
                'pooled_r': meta['pooled_r'],
                'ci_lower_r': meta['ci_lower_r'],
                'ci_upper_r': meta['ci_upper_r'],
                'Q': meta['Q'],
                'I2': meta['I2'],
                'tau2': meta['tau2'],
                'p_heterogeneity': meta['p_heterogeneity'],
                'k': meta['k'],
                'mean_r': np.mean(r_values),  # Simple mean for comparison
            })
            
            # Store forest data for plotting
            for i, (_, row) in enumerate(sub.iterrows()):
                z_i = z_values[i]
                se_i = se_values[i]
                ci_low_r = fisher_z_inv(z_i - 1.96 * se_i)
                ci_high_r = fisher_z_inv(z_i + 1.96 * se_i)
                
                forest_records.append({
                    'pair_name': pair,
                    'category': category,
                    'condition': condition,
                    'dataset': row['dataset'],
                    'r': row[r_col],
                    'n': int(row[n_col]),
                    'ci_lower_r': ci_low_r,
                    'ci_upper_r': ci_high_r,
                    'pooled_r': meta['pooled_r'],
                })
    
    # Save results
    pooled_df = pd.DataFrame(pooled_results)
    pooled_df.to_csv(os.path.join(RESULTS_DIR, 'meta_analysis_pooled.csv'), index=False)
    
    forest_df = pd.DataFrame(forest_records)
    forest_df.to_csv(os.path.join(RESULTS_DIR, 'meta_analysis_forest_data.csv'), index=False)
    
    # Print summary
    print("\n--- Pooled Random-Effects Estimates ---")
    for _, row in pooled_df.iterrows():
        print(f"  {row['pair_name']:20s} {row['condition']:8s}  "
              f"pooled r = {row['pooled_r']:.3f} "
              f"[{row['ci_lower_r']:.3f}, {row['ci_upper_r']:.3f}]  "
              f"I² = {row['I2']:.1f}%  "
              f"(mean r = {row['mean_r']:.3f})")
    
    # Generate forest plots for key ALDH1A1-DA pairs
    key_pairs = ['ALDH1A1-TH', 'ALDH1A1-DDC', 'ALDH1A1-SLC18A2']
    for pair in key_pairs:
        for cond in ['control', 'PD']:
            fdata = forest_df[(forest_df['pair_name'] == pair) & (forest_df['condition'] == cond)]
            if len(fdata) > 0:
                outpath = os.path.join(FIGURES_DIR, f"forest_{pair}_{cond}.pdf")
                make_forest_plot(pair, fdata, cond, outpath)
                print(f"  Saved: {outpath}")
    
    # Z-test comparing control vs PD pooled correlations
    print("\n--- Z-test: Control vs PD Pooled Correlations ---")
    for pair in sorted(target_pairs):
        ctrl_row = pooled_df[(pooled_df['pair_name'] == pair) & (pooled_df['condition'] == 'control')]
        pd_row = pooled_df[(pooled_df['pair_name'] == pair) & (pooled_df['condition'] == 'PD')]
        
        if len(ctrl_row) == 0 or len(pd_row) == 0:
            continue
        
        z_ctrl = fisher_z(ctrl_row['pooled_r'].iloc[0])
        z_pd = fisher_z(pd_row['pooled_r'].iloc[0])
        se_ctrl = ctrl_row.iloc[0]['ci_upper_r']  # Approximate
        se_pd = pd_row.iloc[0]['ci_upper_r']
        
        # Use the pooled SEs directly
        pooled_se_ctrl = (fisher_z(ctrl_row['ci_upper_r'].iloc[0]) - fisher_z(ctrl_row['pooled_r'].iloc[0])) / 1.96
        pooled_se_pd = (fisher_z(pd_row['ci_upper_r'].iloc[0]) - fisher_z(pd_row['pooled_r'].iloc[0])) / 1.96
        
        z_diff = (z_ctrl - z_pd) / np.sqrt(pooled_se_ctrl**2 + pooled_se_pd**2)
        p_diff = 2 * (1 - stats.norm.cdf(abs(z_diff)))
        
        delta_r = pd_row['pooled_r'].iloc[0] - ctrl_row['pooled_r'].iloc[0]
        
        print(f"  {pair:20s}  Δr = {delta_r:+.3f}  z = {z_diff:.2f}  p = {p_diff:.4f}")
    
    print(f"\nResults saved to: {RESULTS_DIR}/meta_analysis_pooled.csv")


if __name__ == '__main__':
    main()
