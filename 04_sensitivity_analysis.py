#!/usr/bin/env python3
"""
04_sensitivity_analysis.py
===========================
Leave-one-out sensitivity analysis and subset analyses for the
ALDH1A1-PD meta-analysis.

Analyses:
  1. Leave-one-out: Recompute mean correlations excluding each dataset
  2. Dataset subset analyses: 4-validated, 5-dataset (excl GSE20292), etc.
  3. SLC6A3 exclusion: ALDH1A1-DA selectivity using only TH, DDC, SLC18A2
  4. 7-dataset vs 6-dataset selectivity comparison

Reads: results/per_dataset_correlations.csv
Produces:
  results/leave_one_out.csv
  results/subset_analyses.csv
"""

import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from config import (
    DATASETS, DECONV_DATASETS, VALIDATED_DATASETS,
    GENE_PAIR_CATEGORIES, get_results_dir
)

RESULTS_DIR = get_results_dir()


def compute_selectivity_from_df(df, datasets=None):
    """Compute selectivity statistics from a correlation DataFrame."""
    if datasets is not None:
        df = df[df['dataset'].isin(datasets)]
    
    aldh1a1 = df[df['category'] == 'ALDH1A1-DA']['delta_r'].dropna()
    dada = df[df['category'] == 'DA-DA']['delta_r'].dropna()
    
    if len(aldh1a1) < 2 or len(dada) < 2:
        return None
    
    mean_a = aldh1a1.mean()
    mean_d = dada.mean()
    sel = mean_a - mean_d
    
    t_stat, t_p = stats.ttest_ind(aldh1a1, dada, equal_var=False)
    u_stat, u_p = stats.mannwhitneyu(aldh1a1, dada, alternative='two-sided')
    
    pooled_std = np.sqrt((aldh1a1.var() * (len(aldh1a1)-1) + dada.var() * (len(dada)-1)) /
                         (len(aldh1a1) + len(dada) - 2))
    d = (mean_a - mean_d) / pooled_std if pooled_std > 0 else np.nan
    
    return {
        'ALDH1A1_mean_dr': mean_a,
        'DADA_mean_dr': mean_d,
        'selectivity': sel,
        't_stat': t_stat, 't_p': t_p,
        'U_stat': u_stat, 'U_p': u_p,
        'cohens_d': d,
        'n_ALDH1A1': len(aldh1a1),
        'n_DADA': len(dada),
    }


def main():
    print("=" * 80)
    print("ALDH1A1-PD META-ANALYSIS: Script 04 — Sensitivity Analysis")
    print("=" * 80)
    
    corr_path = os.path.join(RESULTS_DIR, 'per_dataset_correlations.csv')
    if not os.path.exists(corr_path):
        print("ERROR: per_dataset_correlations.csv not found. Run script 01 first.")
        return
    
    df = pd.read_csv(corr_path)
    target = df[df['category'].isin(['ALDH1A1-DA', 'DA-DA', 'SNCA'])]
    all_datasets = sorted(target['dataset'].unique())
    
    # -----------------------------------------------------------------------
    # 1. Leave-One-Out Analysis
    # -----------------------------------------------------------------------
    print("\n--- Leave-One-Out Sensitivity Analysis ---")
    loo_results = []
    
    for excluded in all_datasets:
        remaining = [d for d in all_datasets if d != excluded]
        sub = target[target['dataset'].isin(remaining)]
        
        for cat in ['ALDH1A1-DA', 'DA-DA', 'SNCA']:
            cat_data = sub[sub['category'] == cat]
            mean_dr = cat_data['delta_r'].mean()
            
            loo_results.append({
                'excluded_dataset': excluded,
                'category': cat,
                'mean_delta_r': mean_dr,
                'n_observations': len(cat_data),
                'n_datasets': len(remaining),
            })
        
        # Selectivity
        sel = compute_selectivity_from_df(sub)
        if sel:
            loo_results.append({
                'excluded_dataset': excluded,
                'category': 'Selectivity',
                'mean_delta_r': sel['selectivity'],
                'n_observations': sel['n_ALDH1A1'] + sel['n_DADA'],
                'n_datasets': len(remaining),
            })
    
    print(f"  Leave-one-out for ALDH1A1-DA mean Δr:")
    for r in loo_results:
        if r['category'] == 'ALDH1A1-DA':
            print(f"    Excl {r['excluded_dataset']}: mean Δr = {r['mean_delta_r']:.3f}")
    
    print(f"\n  Leave-one-out for Selectivity:")
    for r in loo_results:
        if r['category'] == 'Selectivity':
            print(f"    Excl {r['excluded_dataset']}: selectivity = {r['mean_delta_r']:.3f}")
    
    # -----------------------------------------------------------------------
    # 2. Subset Analyses
    # -----------------------------------------------------------------------
    print("\n--- Subset Analyses ---")
    subset_results = []
    
    subsets = {
        '7-dataset (all)': all_datasets,
        '6-dataset (deconv)': DECONV_DATASETS,
        '4-validated': VALIDATED_DATASETS,
        '5-dataset (excl GSE20292)': [d for d in DECONV_DATASETS if d != 'GSE20292'],
    }
    
    for name, ds_list in subsets.items():
        available = [d for d in ds_list if d in all_datasets]
        sel = compute_selectivity_from_df(target, available)
        if sel:
            result = {'subset': name, 'n_datasets': len(available), **sel}
            subset_results.append(result)
            print(f"  {name:30s}: selectivity = {sel['selectivity']:.3f}, "
                  f"t = {sel['t_stat']:.3f}, p(t) = {sel['t_p']:.4f}, "
                  f"U = {sel['U_stat']:.0f}, p(U) = {sel['U_p']:.4f}, "
                  f"d = {sel['cohens_d']:.3f}")
    
    # -----------------------------------------------------------------------
    # 3. SLC6A3 Exclusion Analysis
    # -----------------------------------------------------------------------
    print("\n--- SLC6A3 Exclusion ---")
    
    # Exclude ALDH1A1-SLC6A3 from ALDH1A1-DA pairs
    for name, ds_list in subsets.items():
        available = [d for d in ds_list if d in all_datasets]
        sub = target[target['dataset'].isin(available)]
        
        # ALDH1A1-DA without SLC6A3
        aldh1a1_no_slc6a3 = sub[
            (sub['category'] == 'ALDH1A1-DA') & 
            (~sub['pair_name'].str.contains('SLC6A3'))
        ]['delta_r'].dropna()
        
        if len(aldh1a1_no_slc6a3) > 0:
            mean_dr = aldh1a1_no_slc6a3.mean()
            print(f"  {name}: ALDH1A1-DA (excl SLC6A3) mean Δr = {mean_dr:.3f} "
                  f"(n={len(aldh1a1_no_slc6a3)} observations)")
    
    # -----------------------------------------------------------------------
    # 4. Per-Pair Mean Δr (for Table 2A verification)
    # -----------------------------------------------------------------------
    print("\n--- Per-Pair Mean Δr Across All Datasets ---")
    
    for pair_name in sorted(target['pair_name'].unique()):
        pair_data = target[target['pair_name'] == pair_name]
        mean_ctrl = pair_data['r_ctrl'].mean()
        mean_pd = pair_data['r_pd'].mean()
        mean_dr = pair_data['delta_r'].mean()
        cat = pair_data['category'].iloc[0]
        print(f"  {pair_name:20s} [{cat:12s}] Ctrl: {mean_ctrl:.3f}  PD: {mean_pd:.3f}  Δr: {mean_dr:.3f}")
    
    # -----------------------------------------------------------------------
    # Save Results
    # -----------------------------------------------------------------------
    pd.DataFrame(loo_results).to_csv(os.path.join(RESULTS_DIR, 'leave_one_out.csv'), index=False)
    pd.DataFrame(subset_results).to_csv(os.path.join(RESULTS_DIR, 'subset_analyses.csv'), index=False)
    
    print(f"\nResults saved to {RESULTS_DIR}/")


if __name__ == '__main__':
    main()
