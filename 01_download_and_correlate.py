#!/usr/bin/env python3
"""
01_download_and_correlate.py
============================
Download seven GEO substantia nigra microarray datasets and compute
per-dataset Pearson correlations between target gene pairs in control
and Parkinson's disease groups.

CRITICAL FILTERING REQUIREMENTS
--------------------------------
Two datasets require non-obvious sample filtering to match the
manuscript's reported sample counts (total n=156, 70 ctrl + 86 PD):

  GSE8397 (expected: 15 ctrl + 24 PD = 39)
    - Dual-platform study: each biological sample was profiled on BOTH
      HG-U133A (GPL96) and HG-U133B (GPL97) arrays.
    - Contains multiple brain regions (substantia nigra, frontal cortex, etc.)
    - MUST filter to: (1) GPL96 / A-chip ONLY, AND (2) substantia nigra ONLY
    - Without this filter, you get ~94 rows (duplicated across platforms)

  GSE49036 (expected: 8 ctrl + 15 PD = 23)
    - Contains control, Parkinson's disease, dementia with Lewy bodies (DLB),
      and PD-dementia samples.
    - MUST exclude DLB and PD-dementia; retain ONLY control and pure PD.
    - Without this filter, you get 28+ samples including non-PD pathologies.

All other datasets (GSE7621, GSE20163, GSE20164, GSE20292, GSE20333)
require only standard control-vs-PD group assignment from annotations.

Expected sample counts per dataset:
  GSE7621:   9 ctrl + 16 PD = 25
  GSE8397:  15 ctrl + 24 PD = 39  (A-chip SN only)
  GSE20163:  9 ctrl +  8 PD = 17
  GSE20164:  5 ctrl +  6 PD = 11
  GSE20292: 18 ctrl + 11 PD = 29
  GSE20333:  6 ctrl +  6 PD = 12
  GSE49036:  8 ctrl + 15 PD = 23  (no DLB/dementia)
  TOTAL:    70 ctrl + 86 PD = 156

Output
------
  results/per_dataset_correlations.csv
  results/dataset_summary.csv
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from itertools import combinations

warnings.filterwarnings('ignore')

# Import shared configuration
sys.path.insert(0, os.path.dirname(__file__))
from config import (
    TARGET_GENES, HOUSEKEEPING, DATASETS, GENE_PAIR_CATEGORIES,
    HOUSEKEEPING_PAIRS, get_results_dir, get_category
)

RESULTS_DIR = get_results_dir()


# ---------------------------------------------------------------------------
# GEO Download Helpers
# ---------------------------------------------------------------------------
def download_series_matrix(gse_id):
    """Download and parse a GEO series matrix file.
    
    Returns: (expression_df, sample_annotations)
    """
    import GEOparse
    gse = GEOparse.get_GEO(geo=gse_id, destdir='/tmp/geo_cache', silent=True)
    return gse


def map_probes_to_genes(gse, platform_id):
    """Map probes to gene symbols using platform annotation.
    
    When multiple probes map to the same gene, retains the probe
    with the highest mean expression across all samples.
    """
    gpl = gse.gpls[platform_id]
    
    # Find gene symbol column
    gene_col = None
    for col in ['Gene Symbol', 'GENE_SYMBOL', 'gene_symbol', 'Symbol']:
        if col in gpl.table.columns:
            gene_col = col
            break
    if gene_col is None:
        raise ValueError(f"Cannot find gene symbol column in {platform_id}")
    
    # Build probe-to-gene mapping
    probe_gene = gpl.table.set_index('ID')[gene_col].dropna()
    # Handle multi-gene probes: take first gene
    probe_gene = probe_gene.apply(lambda x: str(x).split('///')[0].strip())
    probe_gene = probe_gene[probe_gene != '']
    probe_gene = probe_gene[probe_gene != '---']
    
    return probe_gene


# ---------------------------------------------------------------------------
# Dataset-Specific Processors
# ---------------------------------------------------------------------------
def process_GSE8397(gse):
    """GSE8397: Dual-platform (GPL96+GPL97), multiple brain regions.
    
    CRITICAL: Must filter to GPL96 (A-chip) AND substantia nigra ONLY.
    Expected: 15 ctrl + 24 PD = 39 samples.
    """
    # Get probe-to-gene mapping for GPL96 only
    probe_gene = map_probes_to_genes(gse, 'GPL96')
    
    samples = []
    for gsm_name, gsm in gse.gsms.items():
        # Filter 1: GPL96 only (A-chip)
        if gsm.metadata.get('platform_id', [''])[0] != 'GPL96':
            continue
        
        # Filter 2: Substantia nigra only
        # Check title and characteristics for brain region
        title = gsm.metadata.get('title', [''])[0].lower()
        chars = ' '.join([str(c) for c in gsm.metadata.get('characteristics_ch1', [])])
        source = ' '.join([str(s) for s in gsm.metadata.get('source_name_ch1', [])]).lower()
        
        region_text = (title + ' ' + chars + ' ' + source).lower()
        
        if 'substantia nigra' not in region_text and 'sn' not in title.split():
            # Also check for common SN abbreviations
            if '_sn_' not in title and not title.endswith('_sn'):
                continue
        
        # Determine group
        group = None
        all_text = (title + ' ' + chars + ' ' + source).lower()
        if 'control' in all_text or 'normal' in all_text:
            group = 'control'
        elif 'parkinson' in all_text or " pd " in all_text or '_pd_' in all_text or all_text.endswith('pd'):
            group = 'PD'
        
        if group is None:
            continue
            
        samples.append({
            'sample_id': gsm_name,
            'group': group,
            'gsm': gsm
        })
    
    # Build expression matrix from GPL96 samples only
    expr_data = {}
    for s in samples:
        table = s['gsm'].table
        if 'VALUE' in table.columns:
            expr_data[s['sample_id']] = table.set_index('ID_REF')['VALUE']
    
    expr_df = pd.DataFrame(expr_data)
    
    # Map probes to genes
    common_probes = expr_df.index.intersection(probe_gene.index)
    expr_df = expr_df.loc[common_probes]
    expr_df.index = probe_gene[common_probes]
    
    # Keep highest-expression probe per gene
    expr_df['mean_expr'] = expr_df.mean(axis=1)
    expr_df = expr_df.sort_values('mean_expr', ascending=False)
    expr_df = expr_df[~expr_df.index.duplicated(keep='first')]
    expr_df = expr_df.drop(columns='mean_expr')
    
    groups = {s['sample_id']: s['group'] for s in samples}
    
    return expr_df, groups


def process_GSE49036(gse):
    """GSE49036: Contains control, PD, DLB, and PD-dementia.
    
    CRITICAL: Must EXCLUDE DLB and PD-dementia samples.
    Retain ONLY control and pure Parkinson's disease.
    Expected: 8 ctrl + 15 PD = 23 samples.
    """
    probe_gene = map_probes_to_genes(gse, 'GPL570')
    
    samples = []
    for gsm_name, gsm in gse.gsms.items():
        title = gsm.metadata.get('title', [''])[0].lower()
        chars = ' '.join([str(c) for c in gsm.metadata.get('characteristics_ch1', [])])
        source = ' '.join([str(s) for s in gsm.metadata.get('source_name_ch1', [])]).lower()
        all_text = (title + ' ' + chars + ' ' + source).lower()
        
        # Exclude DLB and PD-dementia
        if 'lewy bod' in all_text or 'dlb' in all_text:
            continue
        if 'dementia' in all_text:
            continue
        
        # Determine group
        group = None
        if 'control' in all_text or 'normal' in all_text:
            group = 'control'
        elif 'parkinson' in all_text or 'pd' in all_text:
            # Double-check it's not PD-dementia (already excluded above)
            group = 'PD'
        
        if group is None:
            continue
            
        samples.append({
            'sample_id': gsm_name,
            'group': group,
            'gsm': gsm
        })
    
    # Build expression matrix
    expr_data = {}
    for s in samples:
        table = s['gsm'].table
        if 'VALUE' in table.columns:
            expr_data[s['sample_id']] = table.set_index('ID_REF')['VALUE']
    
    expr_df = pd.DataFrame(expr_data)
    
    # Map probes to genes
    common_probes = expr_df.index.intersection(probe_gene.index)
    expr_df = expr_df.loc[common_probes]
    expr_df.index = probe_gene[common_probes]
    
    # Keep highest-expression probe per gene
    expr_df['mean_expr'] = expr_df.mean(axis=1)
    expr_df = expr_df.sort_values('mean_expr', ascending=False)
    expr_df = expr_df[~expr_df.index.duplicated(keep='first')]
    expr_df = expr_df.drop(columns='mean_expr')
    
    groups = {s['sample_id']: s['group'] for s in samples}
    
    return expr_df, groups


def process_standard_dataset(gse, gse_id):
    """Standard processing for GSE7621, GSE20163, GSE20164, GSE20292, GSE20333.
    
    These datasets have straightforward control-vs-PD annotations
    and do not require platform or disease-subtype filtering.
    """
    platform_id = DATASETS[gse_id]['platform']
    probe_gene = map_probes_to_genes(gse, platform_id)
    
    samples = []
    for gsm_name, gsm in gse.gsms.items():
        # Only include samples from the correct platform
        if gsm.metadata.get('platform_id', [''])[0] != platform_id:
            continue
            
        title = gsm.metadata.get('title', [''])[0].lower()
        chars = ' '.join([str(c) for c in gsm.metadata.get('characteristics_ch1', [])])
        source = ' '.join([str(s) for s in gsm.metadata.get('source_name_ch1', [])]).lower()
        all_text = (title + ' ' + chars + ' ' + source).lower()
        
        group = None
        if 'control' in all_text or 'normal' in all_text:
            group = 'control'
        elif 'parkinson' in all_text or " pd" in all_text or '_pd' in all_text:
            group = 'PD'
        
        if group is None:
            continue
            
        samples.append({
            'sample_id': gsm_name,
            'group': group,
            'gsm': gsm
        })
    
    # Build expression matrix
    expr_data = {}
    for s in samples:
        table = s['gsm'].table
        if 'VALUE' in table.columns:
            expr_data[s['sample_id']] = table.set_index('ID_REF')['VALUE']
    
    expr_df = pd.DataFrame(expr_data)
    
    # Map probes to genes
    common_probes = expr_df.index.intersection(probe_gene.index)
    expr_df = expr_df.loc[common_probes]
    expr_df.index = probe_gene[common_probes]
    
    # Keep highest-expression probe per gene
    expr_df['mean_expr'] = expr_df.mean(axis=1)
    expr_df = expr_df.sort_values('mean_expr', ascending=False)
    expr_df = expr_df[~expr_df.index.duplicated(keep='first')]
    expr_df = expr_df.drop(columns='mean_expr')
    
    groups = {s['sample_id']: s['group'] for s in samples}
    
    return expr_df, groups


# ---------------------------------------------------------------------------
# Correlation Computation
# ---------------------------------------------------------------------------
def compute_correlations(expr_df, groups, gene_pairs):
    """Compute Pearson correlations for gene pairs in control and PD groups.
    
    Parameters
    ----------
    expr_df : DataFrame
        Gene × sample expression matrix (gene names as index)
    groups : dict
        sample_id → 'control' or 'PD'
    gene_pairs : list of tuples
        Gene pairs to compute correlations for
    
    Returns
    -------
    list of dicts with columns: gene1, gene2, r_ctrl, r_pd, delta_r,
                                 n_ctrl, n_pd, p_ctrl, p_pd
    """
    ctrl_samples = [s for s, g in groups.items() if g == 'control' and s in expr_df.columns]
    pd_samples = [s for s, g in groups.items() if g == 'PD' and s in expr_df.columns]
    
    results = []
    for g1, g2 in gene_pairs:
        if g1 not in expr_df.index or g2 not in expr_df.index:
            continue
        
        # Control correlation
        x_ctrl = expr_df.loc[g1, ctrl_samples].astype(float).values
        y_ctrl = expr_df.loc[g2, ctrl_samples].astype(float).values
        
        # Remove NaN pairs
        mask_ctrl = ~(np.isnan(x_ctrl) | np.isnan(y_ctrl))
        if mask_ctrl.sum() >= 3:
            r_ctrl, p_ctrl = stats.pearsonr(x_ctrl[mask_ctrl], y_ctrl[mask_ctrl])
            n_ctrl = int(mask_ctrl.sum())
        else:
            r_ctrl, p_ctrl, n_ctrl = np.nan, np.nan, 0
        
        # PD correlation
        x_pd = expr_df.loc[g1, pd_samples].astype(float).values
        y_pd = expr_df.loc[g2, pd_samples].astype(float).values
        
        mask_pd = ~(np.isnan(x_pd) | np.isnan(y_pd))
        if mask_pd.sum() >= 3:
            r_pd, p_pd = stats.pearsonr(x_pd[mask_pd], y_pd[mask_pd])
            n_pd = int(mask_pd.sum())
        else:
            r_pd, p_pd, n_pd = np.nan, np.nan, 0
        
        delta_r = r_pd - r_ctrl if not (np.isnan(r_pd) or np.isnan(r_ctrl)) else np.nan
        
        results.append({
            'gene1': g1,
            'gene2': g2,
            'pair_name': f"{g1}-{g2}",
            'r_ctrl': r_ctrl,
            'r_pd': r_pd,
            'delta_r': delta_r,
            'n_ctrl': n_ctrl,
            'n_pd': n_pd,
            'p_ctrl': p_ctrl,
            'p_pd': p_pd,
        })
    
    return results


# ---------------------------------------------------------------------------
# Main Pipeline
# ---------------------------------------------------------------------------
def main():
    print("=" * 80)
    print("ALDH1A1-PD META-ANALYSIS: Script 01 — Download & Correlate")
    print("=" * 80)
    
    # Define all gene pairs
    target_pairs = list(combinations(TARGET_GENES, 2))
    all_pairs = target_pairs + HOUSEKEEPING_PAIRS
    
    all_results = []
    dataset_summary = []
    
    for gse_id, info in DATASETS.items():
        print(f"\n{'#' * 60}")
        print(f"  Processing {gse_id}...")
        print(f"  Expected: {info['expected_ctrl']} ctrl + {info['expected_pd']} PD")
        print(f"{'#' * 60}")
        
        try:
            gse = download_series_matrix(gse_id)
            
            # Route to appropriate processor
            if gse_id == 'GSE8397':
                expr_df, groups = process_GSE8397(gse)
            elif gse_id == 'GSE49036':
                expr_df, groups = process_GSE49036(gse)
            else:
                expr_df, groups = process_standard_dataset(gse, gse_id)
            
            # Log2 transform if needed (check if values are likely raw)
            max_val = expr_df.max().max()
            if max_val > 100:
                print(f"  Values appear raw (max={max_val:.1f}), applying log2 transform")
                expr_df = np.log2(expr_df.clip(lower=1))
            
            # Count groups
            n_ctrl = sum(1 for g in groups.values() if g == 'control')
            n_pd = sum(1 for g in groups.values() if g == 'PD')
            
            print(f"  Obtained: {n_ctrl} ctrl + {n_pd} PD = {n_ctrl + n_pd} total")
            
            # Validate against expected counts
            if n_ctrl != info['expected_ctrl'] or n_pd != info['expected_pd']:
                print(f"  ⚠️  WARNING: Expected {info['expected_ctrl']} ctrl + "
                      f"{info['expected_pd']} PD, got {n_ctrl} + {n_pd}")
                print(f"  ⚠️  Proceeding, but results may not match manuscript values.")
            else:
                print(f"  ✓ Sample counts match expected values")
            
            # Check target gene availability
            available = [g for g in TARGET_GENES if g in expr_df.index]
            missing = [g for g in TARGET_GENES if g not in expr_df.index]
            print(f"  Target genes available: {len(available)}/{len(TARGET_GENES)}")
            if missing:
                print(f"  Missing: {missing}")
            
            # Compute correlations
            corr_results = compute_correlations(expr_df, groups, all_pairs)
            
            for r in corr_results:
                r['dataset'] = gse_id
                # Assign category
                pair_key = (r['gene1'], r['gene2'])
                pair_key_rev = (r['gene2'], r['gene1'])
                if pair_key in GENE_PAIR_CATEGORIES:
                    r['category'] = GENE_PAIR_CATEGORIES[pair_key]
                elif pair_key_rev in GENE_PAIR_CATEGORIES:
                    r['category'] = GENE_PAIR_CATEGORIES[pair_key_rev]
                else:
                    r['category'] = 'Housekeeping'
            
            all_results.extend(corr_results)
            
            dataset_summary.append({
                'dataset': gse_id,
                'platform': info['platform'],
                'n_ctrl': n_ctrl,
                'n_pd': n_pd,
                'n_total': n_ctrl + n_pd,
                'n_genes': len(expr_df),
                'target_genes_found': len(available),
                'counts_match': n_ctrl == info['expected_ctrl'] and n_pd == info['expected_pd'],
            })
            
        except Exception as e:
            print(f"  ERROR processing {gse_id}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Save results
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(os.path.join(RESULTS_DIR, 'per_dataset_correlations.csv'), index=False)
    
    summary_df = pd.DataFrame(dataset_summary)
    summary_df.to_csv(os.path.join(RESULTS_DIR, 'dataset_summary.csv'), index=False)
    
    # Print summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    total_ctrl = summary_df['n_ctrl'].sum()
    total_pd = summary_df['n_pd'].sum()
    print(f"Total samples: {total_ctrl} ctrl + {total_pd} PD = {total_ctrl + total_pd}")
    print(f"Datasets processed: {len(summary_df)}")
    print(f"Datasets with matching counts: {summary_df['counts_match'].sum()}/{len(summary_df)}")
    print(f"\nResults saved to: {RESULTS_DIR}/per_dataset_correlations.csv")
    
    # Print mean correlations by category
    target_results = results_df[results_df['category'].isin(['ALDH1A1-DA', 'DA-DA', 'SNCA'])]
    print("\n--- Mean Δr by Category (across all datasets) ---")
    for cat in ['ALDH1A1-DA', 'DA-DA', 'SNCA']:
        cat_data = target_results[target_results['category'] == cat]
        mean_dr = cat_data['delta_r'].mean()
        print(f"  {cat}: mean Δr = {mean_dr:.3f} (n={len(cat_data)} observations)")


if __name__ == '__main__':
    main()
