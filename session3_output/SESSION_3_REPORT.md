# Session 3 Report: Differential Abundance + Compartment Re-clustering
## KRT8+ Epithelial Bifurcation Switch in Human Lung Injury

### Dataset
- **Input**: hlca_ready.h5ad (200,000 cells x 25,292 genes)
- **Donors**: 253 total (166 Normal, 67 IPF/ILD, 20 COVID-19)
- **Batches**: HLCA_core (99,720), Adams_ILD_2024 (80,451), Melms_COVID_2021 (19,829)
- **Ambient RNA columns**: None found (status unknown for Melms autopsy samples)

---

## 1. MiloR Differential Abundance Analysis (Main Deliverable)

### Method
Pure Python implementation of MiloR-style differential abundance:
- KNN graph on Harmony-corrected PCA (n_neighbors=30)
- 20,000 neighborhoods (10% sampling, median size 39 cells)
- Negative binomial GLM per neighborhood
- BH FDR correction (significance threshold: FDR < 0.1)

### Results Summary

#### Contrast 1: IPF vs Normal
- **11,449 significant neighborhoods** (of 19,259 tested)
- Enriched in IPF (5,757): Proximal Epithelium, PLIN2+ Fibroblasts, Mesothelial, Myofibroblasts, Proliferating Immune, Endothelial (Macro), Myeloid
- Depleted in IPF (5,692): Airway Epithelium, Alveolar Epithelium, AT2, AT1, Submucosal Gland, B cells, Plasma cells, CD4+ T, Endothelial, Lymphoid

#### Contrast 2: COVID vs Normal
- **2,044 significant neighborhoods** (of 18,195 tested)
- Nearly all enriched in COVID (2,039 up, 5 down)
- Top enriched: Cycling NK/T (logFC 5.5), Macrophages (5.3), Plasma cells (5.2), Monocytes (5.1), CD8+ T (4.8), AT2 (4.6), Fibroblasts (4.5)
- Strong inflammatory cell influx signature

#### Contrast 3: COVID vs IPF (Direct)
- **1,411 significant neighborhoods** (of 12,894 tested)
- Enriched in COVID vs IPF (1,323): Plasma cells (6.3), AT2 (5.3), Macrophages (5.2), CD8+ T (5.4), Cycling NK/T (5.0), Fibroblasts (4.9)
- Depleted in COVID vs IPF (88): Proximal Epithelium, Endothelial (Macro), Mesothelial, PLIN2+ Fibroblasts

### Key Biological Findings
1. **IPF shows massive epithelial remodeling**: Proximal Epithelium expansion with Alveolar/AT2/AT1 loss
2. **COVID shows immune activation**: Macrophage, Plasma cell, CD8+ T influx dominates
3. **Fibrotic niche in IPF**: PLIN2+ Fibroblasts, Mesothelial, Myofibroblasts specifically enriched
4. **Confirms Session 1 descriptive findings** with formal statistical testing

---

## 2. Epithelial Compartment Re-clustering

### Subset
- **86,959 epithelial cells** (9 cell types)
- Normal: 59,414 | IPF/ILD: 25,049 | COVID-19: 2,386
- Includes: Airway Epithelium (36,566), Proximal Epithelium (23,806), Alveolar Epithelium (12,013), Distal Epithelium (9,054), AT2 (1,928), AT1 (1,602), Submucosal Gland (789), Proliferating Epi (788), Other Epi (303)

### Processing
- 3,000 HVGs (batch-aware)
- 50 PCs, Harmony batch correction (converged 2 iterations)
- UMAP (min_dist=0.3), Leiden clustering at res 0.5 (22 clusters) and 1.0 (36 clusters)

### Notable Clusters
- **Cluster 10**: Highly enriched in COVID (1,469 of 3,462 cells = 42%) — likely contains DATP/stress states
- **Cluster 8**: Mixed COVID-enriched (501/3,651 = 14%)
- **Cluster 3**: IPF-dominant (4,579/5,746 = 80%) — likely aberrant basaloid
- **Cluster 0**: IPF-dominant (6,028/7,419 = 81%)
- **Cluster 35**: Small IPF-enriched cluster (49/63)

### Output
- `hlca_epithelial_subset.h5ad` — ready for CellRank (Session 4)
- Marker dotplot: all 24 markers available (KRT8, SFTPC, AGER, KRT17, etc.)

---

## 3. Myeloid Compartment Re-clustering

### Subset
- **62,527 myeloid cells** (5 cell types)
- Normal: 35,962 | IPF/ILD: 22,057 | COVID-19: 4,508
- Includes: Myeloid (57,039), Macrophages (4,436), Monocytes (601), Mast cells (267), DCs (184)

### Processing
- 3,000 HVGs, 50 PCs, Harmony, UMAP, Leiden res 1.0 (26 clusters)

### SPP1+ Macrophage Identification
- **Cluster 0** is the SPP1+ macrophage cluster: mean SPP1 = 1.45 (highest), 7,654 cells
  - IPF-enriched: 6,084/7,654 (79.5%) from IPF, only 549 from COVID
  - This confirms the profibrotic SPP1+ macrophage niche in IPF
- **SPP1 expression by disease**: IPF (0.64) > COVID (0.39) > Normal (0.14)
- Secondary SPP1+ clusters: 19 (mixed), 5 (IPF-enriched)

### Output
- `hlca_myeloid_subset.h5ad` — ready for LIANA+ (Session 7)
- 19/20 myeloid markers available (TNF missing from gene set)

---

## 4. Figures Generated

| Figure | Format | Description |
|--------|--------|-------------|
| milo_beeswarm_ipf_vs_normal | png/svg | Neighborhoods on UMAP, colored by logFC |
| milo_beeswarm_covid_vs_normal | png/svg | Neighborhoods on UMAP, colored by logFC |
| milo_beeswarm_covid_vs_ipf | png/svg | Neighborhoods on UMAP, colored by logFC |
| milo_volcano_all_contrasts | png/svg | logFC vs -log10(FDR), 3 panels |
| milo_celltype_heatmap | png/svg | Mean logFC by cell type across contrasts |
| epi_umap_combined | png/svg | 3-panel: disease, cell type, Leiden |
| epi_umap_disease | png/svg | Disease coloring |
| epi_umap_leiden | png/svg | Leiden 1.0 clusters |
| epi_dotplot_markers | png/svg | 24 KRT8+/DATP/basaloid markers |
| myeloid_umap_combined | png/svg | 3-panel: disease, cell type, Leiden |
| myeloid_umap_disease | png/svg | Disease coloring |
| myeloid_dotplot_markers | png/svg | 19 SPP1/inflammatory markers |

## 5. Tables Generated

| Table | Description |
|-------|-------------|
| milo_da_results_all_contrasts.csv | All 60,000 neighborhood DA results |
| milo_celltype_summary.csv | Cell type-level DA summary |
| milo_nhood_info.csv | Neighborhood UMAP coords + logFC/FDR |
| epi_cluster_composition.csv | Epithelial cluster x disease |
| epi_cluster_celltype.csv | Epithelial cluster x original cell type |
| myeloid_cluster_composition.csv | Myeloid cluster x disease |

## 6. Data Files for Downstream Sessions

| File | Size | Use |
|------|------|-----|
| hlca_epithelial_subset.h5ad | ~87K cells | Session 4: CellRank trajectory |
| hlca_myeloid_subset.h5ad | ~63K cells | Session 7: LIANA+ communication |

---

## 7. Flags for Session 4

1. **Ambient RNA status UNKNOWN** for Melms COVID autopsy samples — no QC columns found
2. **COVID epithelial cells are limited** (2,386) — CellRank may need careful subsetting
3. **Epithelial Cluster 10** (COVID-enriched, 42%) is the primary candidate for DATP/stress state
4. **Gene names**: var_names must be set from `feature_name` column before marker analysis
5. **Harmony embedding**: stored as `X_pca_harmony` in epithelial/myeloid subsets
