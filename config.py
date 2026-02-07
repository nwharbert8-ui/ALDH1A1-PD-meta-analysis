"""
config.py
=========
Shared configuration, constants, and utility functions for the
ALDH1A1-PD meta-analysis pipeline.
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from itertools import combinations

# ---------------------------------------------------------------------------
# Target Genes and Gene Pairs
# ---------------------------------------------------------------------------
TARGET_GENES = ['ALDH1A1', 'TH', 'DDC', 'SLC18A2', 'SLC6A3', 'SNCA']
HOUSEKEEPING = ['GAPDH', 'ACTB', 'RPL13A', 'RPS18', 'HPRT1', 'B2M', 'UBC', 'PPIA']

TARGET_PAIRS = list(combinations(TARGET_GENES, 2))

GENE_PAIR_CATEGORIES = {
    ('ALDH1A1', 'TH'): 'ALDH1A1-DA',
    ('ALDH1A1', 'DDC'): 'ALDH1A1-DA',
    ('ALDH1A1', 'SLC18A2'): 'ALDH1A1-DA',
    ('ALDH1A1', 'SLC6A3'): 'ALDH1A1-DA',
    ('ALDH1A1', 'SNCA'): 'SNCA',
    ('TH', 'DDC'): 'DA-DA',
    ('TH', 'SLC18A2'): 'DA-DA',
    ('TH', 'SLC6A3'): 'DA-DA',
    ('DDC', 'SLC18A2'): 'DA-DA',
    ('DDC', 'SLC6A3'): 'DA-DA',
    ('SLC18A2', 'SLC6A3'): 'DA-DA',
    ('TH', 'SNCA'): 'SNCA',
    ('DDC', 'SNCA'): 'SNCA',
    ('SLC18A2', 'SNCA'): 'SNCA',
    ('SLC6A3', 'SNCA'): 'SNCA',
}

HOUSEKEEPING_PAIRS = [
    ('GAPDH', 'ACTB'),
    ('RPL13A', 'RPS18'),
    ('HPRT1', 'B2M'),
    ('UBC', 'PPIA'),
]

# ---------------------------------------------------------------------------
# Dataset Configuration
# ---------------------------------------------------------------------------
DATASETS = {
    'GSE7621':  {'platform': 'GPL570',  'expected_ctrl': 9,  'expected_pd': 16},
    'GSE8397':  {'platform': 'GPL96',   'expected_ctrl': 15, 'expected_pd': 24},
    'GSE20163': {'platform': 'GPL96',   'expected_ctrl': 9,  'expected_pd': 8},
    'GSE20164': {'platform': 'GPL96',   'expected_ctrl': 5,  'expected_pd': 6},
    'GSE20292': {'platform': 'GPL96',   'expected_ctrl': 18, 'expected_pd': 11},
    'GSE20333': {'platform': 'GPL201',  'expected_ctrl': 6,  'expected_pd': 6},
    'GSE49036': {'platform': 'GPL570',  'expected_ctrl': 8,  'expected_pd': 15},
}

# Datasets included in deconvolution (GSE7621 excluded)
DECONV_DATASETS = ['GSE8397', 'GSE20163', 'GSE20164', 'GSE20292', 'GSE20333', 'GSE49036']

# Datasets with validated DA neuron depletion
VALIDATED_DATASETS = ['GSE8397', 'GSE20163', 'GSE20164', 'GSE49036']

# ---------------------------------------------------------------------------
# Cell Type Signatures (Kamath et al. 2022, Nature Neuroscience)
# Target genes EXCLUDED from ALL signatures (zero circularity)
# ---------------------------------------------------------------------------
CELL_TYPE_SIGNATURES = {
    'DA_Vulnerable_SOX6': [
        'SOX6', 'AGTR1', 'PART1', 'GRM8', 'FAT3',
        'KCNJ6', 'CHRNA6', 'ANKRD20A1', 'RIT2', 'GRIK1',
        'ROBO2', 'PTPRT', 'CDH12', 'DCC', 'TENM1',
        'KIAA1217', 'LSAMP', 'NRG1', 'PRKG1', 'GRID2'
    ],
    'DA_Resistant_CALB1': [
        'CALB1', 'TRHR', 'CALCR', 'LYPD6B', 'PPP1R17',
        'CBLN1', 'GRP', 'SYNPR', 'KCNIP4', 'OTX2',
        'FOXP2', 'EBF1', 'ONECUT2', 'LMO3', 'SOX5'
    ],
    'GABAergic': [
        'GAD1', 'GAD2', 'SLC32A1', 'PVALB', 'SST',
        'VIP', 'LAMP5', 'ADARB2', 'ERBB4', 'LHX6',
        'GRIP1', 'NXPH1', 'SOX14'
    ],
    'Glutamatergic': [
        'SLC17A6', 'SLC17A7', 'SATB2', 'TBR1', 'FEZF2',
        'CUX2', 'RORB', 'THEMIS', 'SLC17A8'
    ],
    'Astrocytes': [
        'AQP4', 'GFAP', 'GJA1', 'SLC1A2', 'SLC1A3',
        'ALDH1L1', 'SOX9', 'NDRG2', 'FABP7', 'GPC5',
        'AGT', 'FGFR3', 'CLU'
    ],
    'Microglia': [
        'CX3CR1', 'P2RY12', 'TMEM119', 'CSF1R', 'HEXB',
        'TREM2', 'C1QA', 'C1QB', 'C1QC', 'ITGAM',
        'AIF1', 'CD68', 'TYROBP'
    ],
    'Oligodendrocytes': [
        'MBP', 'PLP1', 'MOG', 'MAG', 'OPALIN',
        'CNP', 'CLDN11', 'TF', 'MOBP', 'ERMN',
        'ST18', 'FA2H'
    ],
    'Endothelial': [
        'CLDN5', 'FLT1', 'PECAM1', 'VWF', 'ERG',
        'ESAM', 'TIE1', 'CDH5', 'EMCN', 'PODXL'
    ],
}


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
def get_project_root():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def get_results_dir():
    d = os.path.join(get_project_root(), 'results')
    os.makedirs(d, exist_ok=True)
    return d

def get_figures_dir():
    d = os.path.join(get_project_root(), 'figures')
    os.makedirs(d, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Category Lookup Helper
# ---------------------------------------------------------------------------
def get_category(g1, g2):
    """Look up gene pair category from either order."""
    return GENE_PAIR_CATEGORIES.get((g1, g2),
           GENE_PAIR_CATEGORIES.get((g2, g1), 'Other'))
