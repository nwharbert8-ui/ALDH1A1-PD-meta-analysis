# ALDH1A1-PD Meta-Analysis

**ALDH1A1-Dopaminergic Gene Co-expression in Human Substantia Nigra: Meta-Analysis of Disease-Associated Correlation Changes Across Seven Independent Parkinson's Disease Datasets**

Drake Harbert — Inner Architecture LLC

---

## Overview

This repository contains the complete analysis code for a systematic meta-analysis of ALDH1A1-dopaminergic gene co-expression in human substantia nigra across seven independent Parkinson's disease microarray datasets from the Gene Expression Omnibus.

### Key Findings

- **Strong baseline co-expression:** ALDH1A1 shows consistently strong correlations with dopaminergic pathway genes in healthy substantia nigra (mean r > 0.92 for TH, DDC, and SLC18A2)
- **Selective attenuation in PD:** ALDH1A1-dopamine correlations are more attenuated in PD (mean Δr = −0.336) than dopamine-dopamine correlations (mean Δr = −0.143)
- **Selectivity survives cell type adjustment:** After reference-based cell type enrichment analysis, selectivity increases from −0.190 (raw) to −0.210 (adjusted); permutation p = 0.0052
- **Cross-dataset consistency:** The attenuation pattern is observed across all seven independent datasets

## Datasets

| Dataset  | Platform | Ctrl (n) | PD (n) | Total | Notes |
|----------|----------|----------|--------|-------|-------|
| GSE7621  | GPL570   | 9        | 16     | 25    | |
| GSE8397  | GPL96    | 15       | 24     | 39    | A-chip only; SN subset |
| GSE20163 | GPL96    | 9        | 8      | 17    | |
| GSE20164 | GPL96    | 5        | 6      | 11    | |
| GSE20292 | GPL96    | 18       | 11     | 29    | |
| GSE20333 | GPL201   | 6        | 6      | 12    | |
| GSE49036 | GPL570   | 8        | 15     | 23    | PD only; no DLB/dementia |
| **Total**|          | **70**   | **86** | **156** | |

### Critical Filtering Notes

**GSE8397:** This dual-platform study profiled each biological sample on both HG-U133A (GPL96) and HG-U133B (GPL97) arrays, and includes samples from multiple brain regions. The code filters to **GPL96 (A-chip) only** and **substantia nigra only**, yielding 15 ctrl + 24 PD = 39 samples. Without this filter, you would get ~94 rows with duplicated profiling and non-SN tissue.

**GSE49036:** This dataset includes control, Parkinson's disease, dementia with Lewy bodies (DLB), and PD-dementia samples. The code **excludes DLB and PD-dementia**, retaining only control and pure PD samples (8 ctrl + 15 PD = 23). Without this filter, you would get 28+ samples including non-PD pathologies.

## Repository Structure

```
ALDH1A1-PD-meta-analysis/
├── README.md
├── LICENSE                       # MIT License
├── requirements.txt
├── scripts/
│   ├── config.py                 # Shared constants, gene lists, signatures
│   ├── 01_download_and_correlate.py   # Download GEO data + per-dataset correlations
│   ├── 02_meta_analysis.py            # Random-effects meta-analysis (DerSimonian-Laird)
│   ├── 03_cell_type_enrichment.py     # Kamath et al. 2022 deconvolution + permutation test
│   └── 04_sensitivity_analysis.py     # Leave-one-out + subset analyses
├── results/                      # Generated CSV outputs (created by scripts)
└── figures/                      # Generated plots (created by scripts)
```

## Installation

```bash
# Clone repository
git clone https://github.com/YOUR_USERNAME/ALDH1A1-PD-meta-analysis.git
cd ALDH1A1-PD-meta-analysis

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or: venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt
```

### Dependencies

- Python ≥ 3.9
- GEOparse (GEO dataset download and parsing)
- pandas, numpy, scipy (data analysis and statistics)
- matplotlib (figure generation)

## Usage

Run scripts sequentially from the repository root:

```bash
# Step 1: Download datasets and compute per-dataset correlations
python scripts/01_download_and_correlate.py

# Step 2: Random-effects meta-analysis with forest plots
python scripts/02_meta_analysis.py

# Step 3: Cell type enrichment analysis and permutation testing
# (Note: re-downloads datasets; uses cached files if available)
python scripts/03_cell_type_enrichment.py

# Step 4: Sensitivity analyses (leave-one-out, subsets, SLC6A3 exclusion)
python scripts/04_sensitivity_analysis.py
```

Scripts must be run in order (each depends on outputs from previous scripts). The first run downloads ~500 MB of GEO data; subsequent runs use cached files.

## Methods Summary

### Correlation Analysis
Pearson product-moment correlations computed between all pairs of six target genes (ALDH1A1, TH, DDC, SLC18A2, SLC6A3, SNCA) separately for control and PD groups within each dataset. Correlation change Δr = r_PD − r_control.

### Meta-Analysis
Fisher's z-transformation → DerSimonian-Laird random-effects pooling → back-transformation. Heterogeneity assessed via Cochran's Q and I² index.

### Cell Type Enrichment
Reference-based analysis using marker genes from Kamath et al. (2022) single-nucleus RNA-seq of human substantia nigra. **All six target genes excluded from all signatures** (zero circularity). Eight cell type populations scored: DA_Vulnerable (SOX6+/AGTR1+), DA_Resistant (CALB1+), GABAergic, glutamatergic, astrocytes, microglia, oligodendrocytes, endothelial.

### Selectivity Testing
- **Primary test:** Permutation (n=5,000) shuffling disease labels at sample level within datasets, preserving gene pair dependency structure
- **Supplementary:** Welch's t-test, Mann-Whitney U, Cohen's d

## Key Results Summary

### Selectivity Analysis

| Analysis | ALDH1A1 Δr | DA-DA Δr | Selectivity | p (perm) |
|----------|-----------|---------|-------------|----------|
| 6-dataset raw | −0.338 | −0.147 | −0.190 | 0.0052 |
| 7-dataset raw | −0.336 | −0.143 | −0.193 | — |
| 4-validated raw | −0.173 | +0.058 | −0.231 | — |
| 6-dataset adjusted | −0.304 | −0.094 | −0.210 | — |

### Deconvolution Validation

| Dataset | Cohen's d | p-value | Validated |
|---------|----------|---------|-----------|
| GSE8397 | 2.94 | < 0.001 | ✓ |
| GSE20163 | 1.46 | 0.009 | ✓ |
| GSE20164 | 1.86 | 0.014 | ✓ |
| GSE49036 | 1.35 | 0.006 | ✓ |
| GSE20292 | — | > 0.14 | ✗ |
| GSE20333 | — | > 0.14 | ✗ |

## Cell Type Signatures

Marker genes for each cell type were derived from Kamath et al. (2022), *Nature Neuroscience* 25, 588–595. Complete lists are in `scripts/config.py` and Supplementary Table S4 of the manuscript.

**Anti-circularity design:** ALDH1A1, TH, DDC, SLC18A2, SLC6A3, and SNCA were excluded from all cell type signature gene lists. This prevents any mathematical coupling between the deconvolution estimates and the correlation analysis.

## Reproducibility Notes

- All datasets are publicly available from NCBI GEO
- GEOparse caches downloaded files in `/tmp/geo_cache/`
- Expected sample count validation is built into Script 01 (warnings issued if counts don't match)
- Permutation test uses random seed for exact reproducibility (set `np.random.seed(42)` in Script 03 if needed)

## Citation

If you use this code, please cite:

> Harbert, D. (2025). ALDH1A1-Dopaminergic Gene Co-expression in Human Substantia Nigra: Meta-Analysis of Disease-Associated Correlation Changes Across Seven Independent Parkinson's Disease Datasets. *Frontiers in Neuroscience* [under review].

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contact

Drake Harbert — NWHarbert8@gmail.com  
Inner Architecture LLC, Canton, Ohio
