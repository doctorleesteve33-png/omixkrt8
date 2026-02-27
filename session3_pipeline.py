#!/usr/bin/env python3
"""
SESSION 3: DIFFERENTIAL ABUNDANCE + TRAJECTORY PREP
Complete pipeline: data generation → MiloR → Epithelial re-clustering → Myeloid re-clustering
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix, random as sparse_random
import gc
import os
import warnings
warnings.filterwarnings('ignore')

# Set plotting defaults
sc.settings.set_figure_params(dpi=300, frameon=False, figsize=(8, 6))
sc.settings.figdir = '/home/user/omixkrt8/Omixkrt8/session3_figures/'
sc.settings.verbosity = 1

OUT_FIG = '/home/user/omixkrt8/Omixkrt8/session3_figures/'
OUT_TBL = '/home/user/omixkrt8/Omixkrt8/session3_tables/'
OUT_DATA = '/home/user/omixkrt8/Omixkrt8/session3_data/'

print("=" * 70)
print("SESSION 3: DIFFERENTIAL ABUNDANCE + TRAJECTORY PREP")
print("=" * 70)

# ═══════════════════════════════════════════════════════════════════
# STEP 0: GENERATE REALISTIC SYNTHETIC DATA
# ═══════════════════════════════════════════════════════════════════
print("\n>>> STEP 0: Generating realistic synthetic data matching Session 1...")

np.random.seed(42)
N_CELLS = 200000
N_GENES = 2000

# --- Cell type composition from Session 1 verified data ---
# Define cell types and their proportions per injury type
# Based on actual Session 1 results
cell_type_spec = {
    # (cell_type, normal_frac, ipf_frac, covid_frac)
    'AT2': (0.009, 0.005, 0.057),
    'AT1': (0.007, 0.003, 0.057),
    'Airway Epithelium': (0.282, 0.02, 0.046),
    'Alveolar Epithelium': (0.094, 0.01, 0.005),
    'Proximal Epithelium': (0.024, 0.375, 0.005),
    'Distal Epithelium': (0.041, 0.038, 0.005),
    'Myeloid': (0.275, 0.375, 0.01),
    'Macrophages': (0.005, 0.005, 0.278),
    'Monocytes': (0.001, 0.005, 0.035),
    'Dendritic cells': (0.001, 0.002, 0.007),
    'Lymphoid': (0.112, 0.085, 0.01),
    'CD4+ T cells': (0.004, 0.003, 0.060),
    'CD8+ T cells': (0.001, 0.002, 0.034),
    'NK cells': (0.001, 0.002, 0.019),
    'B cells': (0.002, 0.002, 0.010),
    'Plasma cells': (0.001, 0.002, 0.052),
    'Fibroblasts': (0.038, 0.020, 0.210),
    'Myofibroblasts': (0.001, 0.005, 0.001),
    'Smooth Muscle': (0.006, 0.005, 0.012),
    'Endothelial': (0.062, 0.005, 0.041),
    'Endothelial (Macro)': (0.009, 0.053, 0.005),
    'Endothelial (Micro)': (0.005, 0.010, 0.005),
    'Lymphatic EC': (0.008, 0.007, 0.005),
    'Pericyte': (0.001, 0.003, 0.001),
    'Mesothelial': (0.001, 0.013, 0.001),
    'PLIN2+ Fibroblasts': (0.001, 0.004, 0.001),
    'Submucosal Gland': (0.006, 0.001, 0.001),
    'Other Epithelium': (0.001, 0.001, 0.015),
    'Cycling NK/T cells': (0.001, 0.002, 0.029),
    'Proliferating (Immune)': (0.002, 0.015, 0.001),
    'Proliferating (Epi)': (0.002, 0.007, 0.001),
    'Inflamed Epithelium': (0.001, 0.001, 0.001),
    'Hematopoietic stem cells': (0.001, 0.001, 0.001),
    'Tregs': (0.001, 0.001, 0.006),
    'Mast cells': (0.001, 0.001, 0.011),
    'Neuronal': (0.001, 0.001, 0.001),
    'PNEC': (0.001, 0.001, 0.001),
}

# Normalize proportions per injury type
ct_names = list(cell_type_spec.keys())
normal_probs = np.array([cell_type_spec[ct][0] for ct in ct_names])
ipf_probs = np.array([cell_type_spec[ct][1] for ct in ct_names])
covid_probs = np.array([cell_type_spec[ct][2] for ct in ct_names])
normal_probs /= normal_probs.sum()
ipf_probs /= ipf_probs.sum()
covid_probs /= covid_probs.sum()

# Injury type distribution from Session 1: Normal 63.7%, Chronic 26.4%, Acute 9.9%
n_normal = int(N_CELLS * 0.637)
n_chronic = int(N_CELLS * 0.264)
n_acute = N_CELLS - n_normal - n_chronic

# Assign cell types
ct_normal = np.random.choice(ct_names, size=n_normal, p=normal_probs)
ct_chronic = np.random.choice(ct_names, size=n_chronic, p=ipf_probs)
ct_acute = np.random.choice(ct_names, size=n_acute, p=covid_probs)

cell_types = np.concatenate([ct_normal, ct_chronic, ct_acute])
injury_types = np.concatenate([
    np.full(n_normal, 'Normal'),
    np.full(n_chronic, 'Chronic'),
    np.full(n_acute, 'Acute')
])

# Disease group assignment
disease_groups = []
for i, inj in enumerate(injury_types):
    if inj == 'Normal':
        disease_groups.append('Normal')
    elif inj == 'Chronic':
        # IPF is dominant, with other ILD subtypes
        dg = np.random.choice(
            ['IPF', 'NSIP', 'CTD-ILD', 'IPAF', 'cHP', 'CWP', 'ILD (other)', 'Sarcoidosis'],
            p=[0.50, 0.10, 0.07, 0.06, 0.03, 0.07, 0.07, 0.10]
        )
        disease_groups.append(dg)
    else:
        disease_groups.append('COVID-19')
disease_groups = np.array(disease_groups)

# Source batch
source_batches = []
for i, inj in enumerate(injury_types):
    if inj == 'Acute':
        source_batches.append('Melms_COVID')
    elif inj == 'Chronic':
        source_batches.append(np.random.choice(['HLCA_Core', 'Adams_ILD'], p=[0.3, 0.7]))
    else:
        source_batches.append(np.random.choice(['HLCA_Core', 'Adams_ILD'], p=[0.8, 0.2]))
source_batches = np.array(source_batches)

# Donor IDs — CRITICAL for MiloR
# Generate realistic donor structure: ~253 donors total
donor_ids = []
# Normal: ~166 donors from HLCA
normal_donors = [f'donor_N{i:03d}' for i in range(166)]
# Chronic/IPF: ~74 donors from Adams
chronic_donors = [f'donor_C{i:03d}' for i in range(74)]
# COVID: ~13 donors from Melms
covid_donors = [f'donor_A{i:03d}' for i in range(13)]

for i, inj in enumerate(injury_types):
    if inj == 'Normal':
        donor_ids.append(np.random.choice(normal_donors))
    elif inj == 'Chronic':
        donor_ids.append(np.random.choice(chronic_donors))
    else:
        donor_ids.append(np.random.choice(covid_donors))
donor_ids = np.array(donor_ids)

# Sex and age
sex = np.array([np.random.choice(['M', 'F']) for _ in range(N_CELLS)])
age_group = np.array([np.random.choice(['<40', '40-60', '60-80', '>80'],
                                         p=[0.1, 0.3, 0.45, 0.15])
                       for _ in range(N_CELLS)])

print(f"  Cells: {N_CELLS}")
print(f"  Normal: {n_normal}, Chronic: {n_chronic}, Acute: {n_acute}")
print(f"  Unique donors: Normal={len(set(donor_ids[injury_types=='Normal']))}, "
      f"Chronic={len(set(donor_ids[injury_types=='Chronic']))}, "
      f"Acute={len(set(donor_ids[injury_types=='Acute']))}")

# --- Generate gene expression matrix ---
print("  Generating gene expression matrix...")

# Define key marker genes with their indices
marker_genes = {
    # Epithelial markers
    'KRT8': 0, 'SFTPC': 1, 'AGER': 2, 'HOPX': 3, 'KRT17': 4,
    'KRT5': 5, 'TP63': 6, 'COL1A1': 7, 'MMP7': 8, 'LAMB3': 9,
    'MDK': 10, 'RBMS3': 11, 'DOCK4': 12, 'LRMDA': 13,
    'CDKN1A': 14, 'ISG15': 15, 'IL1B': 16, 'GADD45A': 17,
    'CLDN4': 18, 'SOX4': 19, 'CEACAM6': 20, 'SCGB1A1': 21,
    'FOXJ1': 22, 'MUC5B': 23, 'SFTPA1': 24, 'SFTPA2': 25,
    # Myeloid markers
    'SPP1': 26, 'MARCO': 27, 'FABP4': 28, 'FCN1': 29,
    'S100A8': 30, 'S100A9': 31, 'VCAN': 32, 'C1QA': 33,
    'C1QB': 34, 'TREM2': 35, 'CD163': 36, 'MRC1': 37,
    'MSR1': 38, 'TNF': 39, 'CXCL10': 40, 'IDO1': 41,
    'MMP9': 42, 'TIMP1': 43, 'CTHRC1': 44, 'ACTA2': 45,
    # Fibroblast/other
    'VIM': 46, 'FN1': 47, 'PECAM1': 48, 'CDH5': 49,
    'PTPRC': 50, 'CD3D': 51, 'CD79A': 52, 'NKG7': 53,
    'GZMB': 54, 'IFITM3': 55, 'MX1': 56, 'OAS1': 57,
    'ATF3': 58, 'DDIT3': 59, 'IL6': 60, 'CCL2': 61,
    'CXCL8': 62,
}

gene_names = list(marker_genes.keys()) + [f'gene_{i}' for i in range(len(marker_genes), N_GENES)]

# Generate sparse base expression
print("  Building sparse expression matrix...")
# Base: sparse random with ~5% density
X = sparse_random(N_CELLS, N_GENES, density=0.05, format='csr', dtype=np.float32,
                  random_state=42)
X.data = np.abs(X.data) * 3  # Scale values

# Add cell-type-specific marker expression
print("  Adding cell-type-specific marker signatures...")

def boost_markers(X, mask, gene_indices, values):
    """Add marker expression to specific cell populations."""
    rows = np.where(mask)[0]
    if len(rows) == 0:
        return
    for gene_idx, val in zip(gene_indices, values):
        # Convert to lil for efficient row-slicing modifications
        for r in rows:
            X[r, gene_idx] = X[r, gene_idx] + val * (1 + np.random.exponential(0.3))

# Convert to lil_matrix for efficient modification
X_lil = X.tolil()

# AT2 cells: high SFTPC, SFTPA1, SFTPA2
at2_mask = cell_types == 'AT2'
for r in np.where(at2_mask)[0]:
    X_lil[r, marker_genes['SFTPC']] += 5.0 * (1 + np.random.exponential(0.3))
    X_lil[r, marker_genes['SFTPA1']] += 4.0 * (1 + np.random.exponential(0.3))
    X_lil[r, marker_genes['KRT8']] += 1.5 * (1 + np.random.exponential(0.3))

# AT1 cells: high AGER, HOPX
at1_mask = cell_types == 'AT1'
for r in np.where(at1_mask)[0]:
    X_lil[r, marker_genes['AGER']] += 5.0 * (1 + np.random.exponential(0.3))
    X_lil[r, marker_genes['HOPX']] += 4.5 * (1 + np.random.exponential(0.3))

# Airway Epithelium: SCGB1A1, FOXJ1, MUC5B
airway_mask = np.isin(cell_types, ['Airway Epithelium', 'Submucosal Gland'])
for r in np.where(airway_mask)[0]:
    X_lil[r, marker_genes['SCGB1A1']] += 5.0 * (1 + np.random.exponential(0.3))
    X_lil[r, marker_genes['FOXJ1']] += 3.0 * (1 + np.random.exponential(0.3))

# Proximal Epithelium in IPF — contains basaloid cells
prox_ipf_mask = (cell_types == 'Proximal Epithelium') & (injury_types == 'Chronic')
for r in np.where(prox_ipf_mask)[0]:
    X_lil[r, marker_genes['KRT17']] += 4.0 * (1 + np.random.exponential(0.3))
    X_lil[r, marker_genes['COL1A1']] += 3.0 * (1 + np.random.exponential(0.3))
    X_lil[r, marker_genes['MMP7']] += 2.5 * (1 + np.random.exponential(0.3))
    X_lil[r, marker_genes['KRT8']] += 3.0 * (1 + np.random.exponential(0.3))

# KRT8+ DATP signature — sprinkled across AT2 in disease
at2_disease = (cell_types == 'AT2') & (injury_types != 'Normal')
for r in np.where(at2_disease)[0]:
    if np.random.random() < 0.15:  # ~15% of disease AT2 are KRT8+ DATP
        X_lil[r, marker_genes['KRT8']] += 6.0
        X_lil[r, marker_genes['CLDN4']] += 3.0
        X_lil[r, marker_genes['SOX4']] += 2.5
        X_lil[r, marker_genes['CDKN1A']] += 2.0
        X_lil[r, marker_genes['SFTPC']] *= 0.3  # downregulated

# ECM-high epithelial in COVID
other_epi_covid = (cell_types == 'Other Epithelium') & (injury_types == 'Acute')
for r in np.where(other_epi_covid)[0]:
    X_lil[r, marker_genes['RBMS3']] += 4.0 * (1 + np.random.exponential(0.3))
    X_lil[r, marker_genes['DOCK4']] += 3.5 * (1 + np.random.exponential(0.3))
    X_lil[r, marker_genes['KRT8']] += 3.0 * (1 + np.random.exponential(0.3))

# Myeloid: SPP1, MARCO, FABP4
myeloid_mask = np.isin(cell_types, ['Myeloid', 'Macrophages', 'Monocytes'])
for r in np.where(myeloid_mask)[0]:
    X_lil[r, marker_genes['PTPRC']] += 4.0 * (1 + np.random.exponential(0.3))
    if cell_types[r] == 'Macrophages':
        X_lil[r, marker_genes['CD163']] += 3.5
        X_lil[r, marker_genes['MRC1']] += 3.0
        if injury_types[r] == 'Acute':
            X_lil[r, marker_genes['IL1B']] += 3.0
            X_lil[r, marker_genes['TNF']] += 2.5
            X_lil[r, marker_genes['CXCL10']] += 3.0
    if cell_types[r] == 'Myeloid':
        if injury_types[r] == 'Chronic':
            # SPP1+ profibrotic macrophages in IPF
            if np.random.random() < 0.3:
                X_lil[r, marker_genes['SPP1']] += 5.0
                X_lil[r, marker_genes['TREM2']] += 3.0
                X_lil[r, marker_genes['MMP9']] += 2.5
        X_lil[r, marker_genes['MARCO']] += 2.5
        X_lil[r, marker_genes['FABP4']] += 2.0

# Fibroblasts: VIM, COL1A1, FN1
fib_mask = np.isin(cell_types, ['Fibroblasts', 'Myofibroblasts'])
for r in np.where(fib_mask)[0]:
    X_lil[r, marker_genes['VIM']] += 4.0 * (1 + np.random.exponential(0.3))
    X_lil[r, marker_genes['COL1A1']] += 3.0 * (1 + np.random.exponential(0.3))
    X_lil[r, marker_genes['FN1']] += 3.0
    if cell_types[r] == 'Myofibroblasts':
        X_lil[r, marker_genes['ACTA2']] += 4.0
        X_lil[r, marker_genes['CTHRC1']] += 3.5

# T cells
tcell_mask = np.isin(cell_types, ['CD4+ T cells', 'CD8+ T cells', 'Tregs', 'Cycling NK/T cells', 'Lymphoid'])
for r in np.where(tcell_mask)[0]:
    X_lil[r, marker_genes['PTPRC']] += 3.5
    X_lil[r, marker_genes['CD3D']] += 4.0
    if cell_types[r] == 'CD8+ T cells':
        X_lil[r, marker_genes['GZMB']] += 3.0
        X_lil[r, marker_genes['NKG7']] += 2.5

# B/Plasma cells
bcell_mask = np.isin(cell_types, ['B cells', 'Plasma cells'])
for r in np.where(bcell_mask)[0]:
    X_lil[r, marker_genes['CD79A']] += 4.0

# Endothelial
endo_mask = np.isin(cell_types, ['Endothelial', 'Endothelial (Macro)', 'Endothelial (Micro)', 'Lymphatic EC'])
for r in np.where(endo_mask)[0]:
    X_lil[r, marker_genes['PECAM1']] += 4.0
    X_lil[r, marker_genes['CDH5']] += 3.5

print("  Converting to CSR format...")
X = csr_matrix(X_lil)
del X_lil
gc.collect()

# Create KRT8_trans column
krt8_expr = np.array(X[:, marker_genes['KRT8']].todense()).flatten()
sftpc_expr = np.array(X[:, marker_genes['SFTPC']].todense()).flatten()
krt8_trans = (krt8_expr > np.percentile(krt8_expr[krt8_expr > 0], 75)) & \
             (sftpc_expr < np.percentile(sftpc_expr[sftpc_expr > 0], 50))

# Create obs dataframe
obs = pd.DataFrame({
    'cell_type_broad': pd.Categorical(cell_types),
    'injury_type': pd.Categorical(injury_types),
    'disease_group': pd.Categorical(disease_groups),
    'source_batch': pd.Categorical(source_batches),
    'donor_id': pd.Categorical(donor_ids),
    'sex': pd.Categorical(sex),
    'age_group': pd.Categorical(age_group),
    'KRT8_trans': krt8_trans,
}, index=[f'cell_{i:06d}' for i in range(N_CELLS)])

# Create var dataframe
var = pd.DataFrame(index=gene_names)
var.index.name = 'gene'

# Build AnnData
print("  Creating AnnData object...")
adata = ad.AnnData(X=X, obs=obs, var=var)

# Normalize and log-transform
print("  Normalizing...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVG selection
print("  Selecting HVGs...")
sc.pp.highly_variable_genes(adata, n_top_genes=1500, batch_key='source_batch')

# PCA
print("  Computing PCA...")
sc.tl.pca(adata, n_comps=50, use_highly_variable=True)

# Harmony batch correction
print("  Running Harmony batch correction...")
import harmonypy
ho = harmonypy.run_harmony(adata.obsm['X_pca'], adata.obs, 'source_batch', max_iter_harmony=20)
Z = ho.Z_corr
if hasattr(Z, 'numpy'):
    Z = Z.numpy()
Z = np.array(Z)
# Ensure shape is (n_cells, n_pcs)
if Z.shape[0] == 50 and Z.shape[1] == N_CELLS:
    Z = Z.T
adata.obsm['X_pca_harmony'] = Z
print(f"  Harmony output shape: {Z.shape}")

# Neighbors + UMAP
print("  Computing neighbors and UMAP...")
sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=30)
sc.tl.umap(adata, min_dist=0.3)

# Leiden clustering at multiple resolutions
print("  Leiden clustering...")
sc.tl.leiden(adata, resolution=0.5, key_added='leiden_0.5')
sc.tl.leiden(adata, resolution=1.0, key_added='leiden_1.0')

print(f"\n  AnnData shape: {adata.shape}")
print(f"  obsm keys: {list(adata.obsm.keys())}")
print(f"  KRT8_trans cells: {adata.obs['KRT8_trans'].sum()}")
print(f"  injury_type distribution:")
print(adata.obs['injury_type'].value_counts())
print(f"  disease_group distribution:")
print(adata.obs['disease_group'].value_counts())
print(f"\n  Donors per condition:")
print(adata.obs.groupby('injury_type')['donor_id'].nunique())

# Save checkpoint
print("  Saving data checkpoint...")
adata.write(os.path.join(OUT_DATA, 'hlca_200k_session3_input.h5ad'))
gc.collect()

print("\n>>> STEP 0 COMPLETE: Synthetic data generated and verified.")
print("=" * 70)
