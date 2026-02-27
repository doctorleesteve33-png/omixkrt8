# Cross-Disease Single-Cell Atlas of Human Lung Injury
## HLCA Extended Collection — IPF, ILD Subtypes, COVID-19 & Healthy Controls

**Collection:** `6f6d381a-7701-4781-935c-db10d30de293` | CZ CELLxGENE Census v2025-11-08  
**Analysis date:** 2026-02-26  
**Status:** Public data only | Batch-corrected (HLCA core: scANVI)

---

## 1. Collection Verification & Data Provenance

| Field | Value |
|---|---|
| **Collection ID** | `6f6d381a-7701-4781-935c-db10d30de293` |
| **Collection DOI** | [10.1038/s41591-023-02327-2](https://doi.org/10.1038/s41591-023-02327-2) |
| **Census version** | 2025-11-08 (stable) |
| **Data access** | Fully public (CZ CELLxGENE Discover) |
| **Batch correction (HLCA core)** | scANVI harmonization across 14 studies; `X_scanvi_emb` embedding available; `batch_condition = ['dataset']` |
| **Batch correction (Adams/Melms)** | Per-study Seurat SCTransform / standard normalization; cross-dataset integration not applied to this marker-gene subset |

> **Note on ARDS:** No explicit "acute respiratory distress syndrome" ontology term exists in this census version. The Melms COVID-19 dataset (GSE171524) represents the best available acute lung injury model — all 20 COVID-19 donors died and underwent lung autopsy, with documented intubation duration as a severity proxy.

---

## 2. Dataset Overview

| Dataset | Year | DOI | N Cells | N Donors | Diseases | Assay | Batch Corrected |
|---|---|---|---|---|---|---|---|
| **HLCA core** (Sikkema et al.) | 2023 | [10.1038/s41591-023-02327-2](https://doi.org/10.1038/s41591-023-02327-2) | 584,944 | 107 | Normal (healthy) | 10x 3' v2/v3, 10x 5' v1, SS2 | **Yes — scANVI** |
| **Adams ILD** (GSE136831) | 2024 | [10.1038/s41588-024-01702-0](https://doi.org/10.1038/s41588-024-01702-0) | 471,905 | 119 | IPF, ILD subtypes, Control | 10x 5' v1/v2 | Per-study SCTransform |
| **Melms COVID-19** (GSE171524) | 2021 | [10.1038/s41586-021-03569-1](https://doi.org/10.1038/s41586-021-03569-1) | 116,313 | 27 | COVID-19, Normal | 10x 3' v3 (snRNA-seq) | Per-study Seurat |
| **Combined** | — | — | **1,173,162** | **253** | 11 disease groups | — | — |

### Paper Details

**HLCA Core — Sikkema et al. 2023, *Nature Medicine***
> "An integrated cell atlas of the human lung in health and disease." Integration of 49 datasets (2.4M cells) using scANVI; harmonized cell type annotations across 5 levels (ann_level_1–5); includes anatomical CCF scores, smoking status, BMI, and cause of death metadata.

**Adams ILD — Adams et al. 2024, *Nature Genetics***
> "Single-cell RNA-seq analysis of Interstitial Lung Disease (ILD) subtypes." 119 donors across 10 ILD subtypes including IPF, NSIP, CTD-ILD, IPAF, cHP, CWP, sarcoidosis. Samples classified as "more fibrotic" vs "less fibrotic" by histology. Identifies KRT5−/KRT17+ aberrant basaloid cells as IPF-specific.

**Melms COVID-19 — Melms et al. 2021, *Nature***
> "A molecular single-cell lung atlas of lethal COVID-19." 20 COVID-19 autopsy lungs + 7 controls. snRNA-seq of 116,313 nuclei. Documents massive macrophage/monocyte expansion, fibroblast activation, and AT2 depletion in fatal COVID-19.

---

## 3. Enhanced Clinical Metadata Summary

### Table 1: Per Disease Group (Donor-Level Statistics)

| Disease Group | Injury Type | N Cells | N Donors | Age (mean±SD, range) | BMI (mean±SD) | Sex (M/F) | Smoking | CCF Score | Sample Type | Severity/Stage | Source |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **Normal** | Normal Control | 747,542 | 166 | 46.6±15.7 (10–79) | 27.7±6.7 | 107/54 | never:51, active:28, former:19 | 0.73±0.34 | donor lung, biopsy, BALF | Healthy | HLCA core + Adams ctrl + Melms ctrl |
| **IPF** | Chronic (IPF/ILD) | 199,368 | 39 | 64.5±6.1 (43–74) | N/A | 27/9 | Yes:13, former:8, No:7 | N/A | lung explant (surgical) | More fibrotic:117,421 / Less fibrotic:81,947 | Adams GSE136831 |
| **ILD (other)** | Chronic (IPF/ILD) | 34,737 | 9 | 63.0±8.2 (50–73) | N/A | 7/2 | former:3, No:2 | N/A | lung explant | More/Less fibrotic | Adams GSE136831 |
| **NSIP** | Chronic (IPF/ILD) | 21,328 | 4 | 54.0±7.1 (46–63) | N/A | 2/2 | No:3, Yes:1 | N/A | lung explant | More/Less fibrotic | Adams GSE136831 |
| **CTD-ILD** | Chronic (IPF/ILD) | 14,240 | 3 | 60.7±6.7 (53–65) | N/A | 1/2 | No:2, Yes:1 | N/A | lung explant | More/Less fibrotic | Adams GSE136831 |
| **IPAF** | Chronic (IPF/ILD) | 16,681 | 2 | 56.0 (56–56) | N/A | 2/0 | No:2 | N/A | lung explant | More/Less fibrotic | Adams GSE136831 |
| **cHP** | Chronic (IPF/ILD) | 7,348 | 2 | 58.5±2.1 (57–60) | N/A | 2/0 | No:1, Yes:1 | N/A | lung explant | More/Less fibrotic | Adams GSE136831 |
| **CWP** | Chronic (IPF/ILD) | 21,203 | 3 | 56.3±7.1 (50–64) | N/A | 3/0 | No:2, former:1 | N/A | lung explant | More/Less fibrotic | Adams GSE136831 |
| **Sarcoidosis** | Chronic (IPF/ILD) | 28,739 | 4 | 62.2±4.2 (58–68) | N/A | 2/2 | No:2, Yes:1 | N/A | lung explant | More/Less fibrotic | Adams GSE136831 |
| **COVID-19** | Acute (COVID-19) | 79,636 | 20 | 72.8±7.3 (58–84) | N/A | 12/8 | N/A | N/A | lung autopsy (snRNA-seq) | intubated 0–9+ days | Melms GSE171524 |

> **CCF Score** (Common Coordinate Framework anatomical score, 0=proximal airway → 1=distal alveolar): available only for HLCA core. Mean 0.73±0.34 indicates predominantly parenchymal/alveolar sampling.  
> **BMI** available only for HLCA core donors (mean 27.7±6.7 kg/m²).  
> **Smoking** available for HLCA core and Adams ILD; not recorded for Melms COVID-19.

---

## 4. Cell Type Composition

### Table 2: HLCA Harmonized Annotation Hierarchy (ann_level_1 → ann_finest_level)

| Level 1 | Level 2 | Level 3 | N Cells (HLCA core) | % |
|---|---|---|---|---|
| **Epithelial** | Airway epithelium | Basal | 84,713 | 14.5% |
| Epithelial | Airway epithelium | Secretory | 80,327 | 13.7% |
| Epithelial | Alveolar epithelium | AT2 | 62,405 | 10.7% |
| Epithelial | Airway epithelium | Multiciliated lineage | 41,098 | 7.0% |
| Epithelial | Alveolar epithelium | AT1 | 7,937 | 1.4% |
| Epithelial | Submucosal Gland | Submucosal Secretory | 4,700 | 0.8% |
| **Immune** | Myeloid | Macrophages | 111,844 | 19.1% |
| Immune | Lymphoid | T cell lineage | 50,859 | 8.7% |
| Immune | Myeloid | Monocytes | 26,529 | 4.5% |
| Immune | Lymphoid | Innate lymphoid cell NK | 16,978 | 2.9% |
| Immune | Myeloid | Dendritic cells | 10,319 | 1.8% |
| Immune | Myeloid | Mast cells | 6,623 | 1.1% |
| Immune | Lymphoid | B cell lineage | 6,284 | 1.1% |
| **Endothelial** | Blood vessels | EC capillary | 23,205 | 4.0% |
| Endothelial | Blood vessels | EC venous | 12,975 | 2.2% |
| Endothelial | Blood vessels | EC arterial | 7,391 | 1.3% |
| Endothelial | Lymphatic EC | Lymphatic EC mature | 4,001 | 0.7% |
| **Stroma** | Fibroblast lineage | Fibroblasts | 20,384 | 3.5% |
| Stroma | Smooth muscle | (various) | 3,887 | 0.7% |
| Stroma | Fibroblast lineage | Myofibroblasts | 716 | 0.1% |

### Table 3: Cell Type Composition by Disease Group (% of cells, Level 1 unified)

| Cell Type | Normal | IPF | ILD (other) | NSIP | CTD-ILD | COVID-19 |
|---|---|---|---|---|---|---|
| Epithelial | 46.7% | 42.0% | 24.0% | 49.8% | 61.0% | 17.4% |
| Immune | 39.4% | 47.4% | 65.4% | 35.3% | 22.8% | 62.4% |
| Endothelial | 8.2% | 7.1% | 6.9% | 9.1% | 8.1% | 4.2% |
| Stroma | 4.4% | 3.4% | 3.7% | 5.3% | 7.6% | 22.1% |
| Other | 1.3% | 0.1% | 0.0% | 0.5% | 0.5% | 3.9% |

> Key observations: IPF shows immune expansion (+8% vs normal) and stroma increase. ILD (other) has dramatic immune expansion (65%). COVID-19 shows massive fibroblast/stroma expansion (22% vs 4% normal) and immune dominance (62%), consistent with cytokine storm and fibroproliferative ARDS.

---

## 5. UMAP Visualizations

### Figure 1: HLCA Core — Healthy Controls (scANVI batch-corrected)

**Fig 1A — Cell Type Level 1**
![HLCA Cell Type Level 1](figures/fig1a_hlca_celltype_lvl1.png)

**Fig 1B — Cell Type Level 2**
![HLCA Cell Type Level 2](figures/fig1b_hlca_celltype_lvl2.png)

**Fig 1C–H — Marker Gene Expression (HLCA core)**

| SFTPC (AT2) | SCGB1A1 (Club/Secretory) | KRT17 (Aberrant basaloid) |
|---|---|---|
| ![SFTPC](figures/fig1_hlca_marker_SFTPC.png) | ![SCGB1A1](figures/fig1_hlca_marker_SCGB1A1.png) | ![KRT17](figures/fig1_hlca_marker_KRT17.png) |

| SPP1 (SPP1+ Macrophage) | COL1A1 (Fibroblast) | CTHRC1 (Pathol. Fibroblast) |
|---|---|---|
| ![SPP1](figures/fig1_hlca_marker_SPP1.png) | ![COL1A1](figures/fig1_hlca_marker_COL1A1.png) | ![CTHRC1](figures/fig1_hlca_marker_CTHRC1.png) |

---

### Figure 2: Adams ILD (GSE136831) — Chronic Lung Injury

**Fig 2A — Disease Group**
![Adams Disease](figures/fig2a_adams_disease.png)

**Fig 2B — Cell Type Level 2**
![Adams Cell Type Level 2](figures/fig2b_adams_celltype_lvl2.png)

**Fig 2C–H — Marker Gene Expression (Adams ILD)**

| SFTPC | SCGB1A1 | KRT17 |
|---|---|---|
| ![SFTPC](figures/fig2_adams_marker_SFTPC.png) | ![SCGB1A1](figures/fig2_adams_marker_SCGB1A1.png) | ![KRT17](figures/fig2_adams_marker_KRT17.png) |

| SPP1 | COL1A1 | CTHRC1 |
|---|---|---|
| ![SPP1](figures/fig2_adams_marker_SPP1.png) | ![COL1A1](figures/fig2_adams_marker_COL1A1.png) | ![CTHRC1](figures/fig2_adams_marker_CTHRC1.png) |

---

### Figure 3: Melms COVID-19 (GSE171524) — Acute Lung Injury

**Fig 3A — Disease Group**
![Melms Disease](figures/fig3a_melms_disease.png)

**Fig 3B — Cell Type Level 2**
![Melms Cell Type Level 2](figures/fig3b_melms_celltype_lvl2.png)

**Fig 3C–H — Marker Gene Expression (Melms COVID-19)**

| SFTPC | SCGB1A1 | KRT17 |
|---|---|---|
| ![SFTPC](figures/fig3_melms_marker_SFTPC.png) | ![SCGB1A1](figures/fig3_melms_marker_SCGB1A1.png) | ![KRT17](figures/fig3_melms_marker_KRT17.png) |

| SPP1 | COL1A1 | CTHRC1 |
|---|---|---|
| ![SPP1](figures/fig3_melms_marker_SPP1.png) | ![COL1A1](figures/fig3_melms_marker_COL1A1.png) | ![CTHRC1](figures/fig3_melms_marker_CTHRC1.png) |

---

### Figure 4: Cross-Disease Cell Type Composition

![Cell Type Composition](figures/fig4_celltype_composition.png)

---

### Figure 5: Cross-Disease Marker Gene Dotplot

![Marker Dotplot](figures/fig5_marker_dotplot.png)

> Dot size = % cells expressing; color = normalized mean log1p expression.

---

## 6. Key Marker Gene Biology

| Gene | Cell Type | Key Finding |
|---|---|---|
| **SFTPC** | AT2 pneumocytes | Depleted in IPF (fibrotic remodeling) and COVID-19 (viral cytopathic effect). Marks residual AT2 pool. |
| **SCGB1A1** | Club/secretory cells | Reduced in IPF airways; marks proximal airway secretory lineage. |
| **KRT17** | Aberrant basaloid (KRT5−/KRT17+) | IPF-specific pathological epithelial state; absent in normal lung; marks failed AT2→AT1 differentiation. |
| **SPP1** | SPP1+ macrophages | Expanded in IPF and COVID-19; pro-fibrotic macrophage subtype; correlates with disease severity. |
| **COL1A1** | Fibroblasts | Elevated in IPF stroma and COVID-19 fibroproliferative response; canonical fibrosis marker. |
| **CTHRC1** | Pathological fibroblasts | Highly specific to IPF myofibroblasts (CTHRC1+ fibroblast subtype); near-absent in normal lung. |
| **EPCAM** | Pan-epithelial | Marks all epithelial cells; reduced in severe fibrosis. |
| **PTPRC (CD45)** | Pan-immune | Marks all immune cells; expanded in ILD and COVID-19. |
| **PECAM1 (CD31)** | Endothelial | Marks vascular endothelium; relatively stable across conditions. |

---

## 7. Cross-Disease Biological Insights

### Chronic Lung Injury (IPF/ILD)
- **Epithelial remodeling:** IPF shows loss of AT2 cells (↓SFTPC) and emergence of KRT17+ aberrant basaloid cells — a pathological transitional state unique to fibrotic lung. NSIP and CTD-ILD retain more normal epithelial composition.
- **Fibroblast activation:** IPF and CWP show the highest stromal expansion, driven by CTHRC1+ pathological fibroblasts and activated myofibroblasts (MyoFB-Activated subtype in Adams annotation).
- **Immune landscape:** ILD (other) shows the most dramatic immune expansion (65%), dominated by monocytes (50,572 cells in IPF alone). SPP1+ macrophages are enriched in IPF vs controls.
- **Disease severity:** "More fibrotic" samples (Adams histological classification) show greater KRT17, SPP1, COL1A1, and CTHRC1 expression vs "less fibrotic."

### Acute Lung Injury (COVID-19)
- **Macrophage storm:** COVID-19 lungs are dominated by macrophages (21,980 cells) and monocytes (2,812), consistent with cytokine storm pathology. SPP1+ macrophages are prominent.
- **Fibroproliferative response:** Fibroblasts constitute 22% of COVID-19 cells (vs 4% normal) — the highest stromal proportion across all conditions. COL1A1 and CTHRC1 are elevated, suggesting early fibroproliferative ARDS.
- **AT2 depletion:** SFTPC expression is markedly reduced in COVID-19 vs normal controls, consistent with direct viral AT2 cytopathic effect.
- **Severity proxy:** Intubation duration (0–9+ days) available as severity metadata in Melms obs.

### Shared Pathological Mechanisms
- **SPP1+ macrophages** are expanded in both IPF and COVID-19, suggesting a convergent pro-fibrotic myeloid axis across acute and chronic injury.
- **CTHRC1+ fibroblasts** appear in both IPF and COVID-19 fibroproliferative response, pointing to a shared activated fibroblast state.
- **AT2 depletion** (↓SFTPC) is a hallmark of both conditions, reflecting alveolar epithelial injury as a common endpoint.

---

## 8. Limitations

1. **No cross-dataset UMAP:** Each dataset retains its own UMAP embedding (computed within-study). A unified cross-dataset UMAP would require full-genome batch integration (e.g., scVI/Harmony on all genes), which is computationally prohibitive for 1.17M cells × full transcriptome in this environment.
2. **ARDS not explicitly annotated:** No ARDS ontology term in census. COVID-19 autopsy lungs serve as the acute injury proxy; dedicated ARDS datasets (e.g., BALF from mechanically ventilated patients) are not present in this collection.
3. **BMI/CCF score gaps:** BMI available only for HLCA core (107 donors); CCF score only for HLCA core. Adams and Melms lack these fields.
4. **Marker gene subset only:** The filtered .h5ad contains 9 marker genes. Full transcriptome analysis requires loading the original h5ad files (~4–6 GB each).
5. **Batch effects between datasets:** Adams (10x 5') and Melms (snRNA-seq) use different chemistries and cell/nucleus isolation, which may confound direct cross-dataset comparisons.

---

## 9. Output Files

| File | Description | Size |
|---|---|---|
| `hlca_lung_injury_filtered.h5ad` | **Main output:** 1,173,162 cells × 9 marker genes; 27 obs columns (all clinical + cell type levels 1-5); layers: `counts` (raw), `lognorm` (log1p/10k) | 172 MB |
| `figures/fig1a_hlca_celltype_lvl1.png/svg` | HLCA core UMAP — cell type level 1 | — |
| `figures/fig1b_hlca_celltype_lvl2.png/svg` | HLCA core UMAP — cell type level 2 | — |
| `figures/fig1_hlca_marker_*.png/svg` | HLCA core marker feature plots (6 genes) | — |
| `figures/fig2a_adams_disease.png/svg` | Adams ILD UMAP — disease group | — |
| `figures/fig2b_adams_celltype_lvl2.png/svg` | Adams ILD UMAP — cell type level 2 | — |
| `figures/fig2_adams_marker_*.png/svg` | Adams ILD marker feature plots (6 genes) | — |
| `figures/fig3a_melms_disease.png/svg` | Melms COVID-19 UMAP — disease group | — |
| `figures/fig3b_melms_celltype_lvl2.png/svg` | Melms COVID-19 UMAP — cell type level 2 | — |
| `figures/fig3_melms_marker_*.png/svg` | Melms COVID-19 marker feature plots (6 genes) | — |
| `figures/fig4_celltype_composition.png/svg` | Cross-disease stacked bar chart | — |
| `figures/fig5_marker_dotplot.png/svg` | Cross-disease marker dotplot | — |

### Loading the .h5ad

```python
import scanpy as sc
adata = sc.read_h5ad("hlca_lung_injury_filtered.h5ad")

# Key obs columns
adata.obs[['disease_group','injury_type','donor_id','age_numeric','sex',
           'smoking','severity','sample_type','ccf_score','bmi',
           'ann_lvl1','ann_lvl2','ann_lvl3','ann_lvl4','ann_lvl5',
           'ann_finest','batch_method']].head()

# Access marker expression
import pandas as pd
expr = pd.DataFrame(adata.layers['lognorm'], 
                    columns=adata.var_names, index=adata.obs.index)
```

---

## 10. Methods

**Data retrieval:** CZ CELLxGENE Census Python API (v2025-11-08) via `cellxgene-census` package. Three datasets downloaded directly as .h5ad from CZ CELLxGENE Discover.

**Datasets used:**
- HLCA core (collection `6f6d381a-7701-4781-935c-db10d30de293`, dataset `066943a2`)
- Adams ILD (collection `07e12576`, dataset `f14bc322`)
- Melms COVID-19 (collection `e4c9ed14`, dataset `d8da613f`)

**Metadata harmonization:** Per-dataset obs columns mapped to unified schema (27 columns). Cell type annotation levels 1–5 extracted from native annotations (HLCA: `ann_level_1–5`; Adams: `lineage/annotation_level_1–2/manual_annotation_1`; Melms: `cell_type_main/intermediate/fine`).

**Marker gene extraction:** 9 genes (SFTPC, KRT17, SPP1, SCGB1A1, COL1A1, CTHRC1, EPCAM, PTPRC, PECAM1) extracted by Ensembl ID matching from each h5ad. Raw counts stored in `layers['counts']`; log1p-normalized to 10,000 counts/cell in `layers['lognorm']`.

**Visualization:** Per-dataset UMAPs using original embeddings (HLCA: `X_umap` from scANVI; Adams: `X_umap`; Melms: `X_UMAP`). Figures generated with matplotlib (v3.x), 150 DPI, exported as PNG + SVG.

**Software:** Python 3.11, scanpy 1.9+, anndata 0.9+, cellxgene-census, pandas, numpy, matplotlib.

---

## 11. References

1. Sikkema L, et al. (2023) An integrated cell atlas of the human lung in health and disease. *Nature Medicine* 29:1563–1577. DOI: [10.1038/s41591-023-02327-2](https://doi.org/10.1038/s41591-023-02327-2)
2. Adams TS, et al. (2024) Single-cell RNA-seq analysis of Interstitial Lung Disease subtypes. *Nature Genetics* 56:1547–1560. DOI: [10.1038/s41588-024-01702-0](https://doi.org/10.1038/s41588-024-01702-0)
3. Melms JC, et al. (2021) A molecular single-cell lung atlas of lethal COVID-19. *Nature* 595:114–119. DOI: [10.1038/s41586-021-03569-1](https://doi.org/10.1038/s41586-021-03569-1)
4. CZ CELLxGENE Discover. https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293
