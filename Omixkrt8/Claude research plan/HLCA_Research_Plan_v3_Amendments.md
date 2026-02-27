# Research Plan v3: Targeted Amendments to v2
## What Changed and Why — With Full Fact-Checking

**This document patches v2. The v2 structure, figures, and session plan remain intact.**
**Read v2 first, then apply these 5 amendments.**

---

## AMENDMENT 1: Split Spatial Pipeline — Visium vs. Xenium
### Status: ✅ ADOPTED (Gemini is correct, this was a real error in v2)

**The problem in v2:** Phase 7 (Section 3.7.1) proposed running cell2location on both Visium and Xenium data. This is a methodological error. Visium captures 50 µm spots containing 5–15 cells each, requiring deconvolution to estimate cell type proportions per spot. Xenium is subcellular-resolution imaging-based spatial transcriptomics — the cells are already individually segmented. Running deconvolution on pre-segmented single cells is technically meaningless and reviewers will flag it immediately.

**The v3 fix:** Split the spatial pipeline into two distinct tracks.

### Track A: Visium Data (Deconvolution-Based)
**Datasets**: Franzén 2024 (E-MTAB-14121), Mothes 2023 (Zenodo 7533977)
**Method**: cell2location (Kleshchevnikov 2022)
**Purpose**: Map broad regional gradients — show how cell type composition transitions from normal parenchyma → transitional zone → fibrotic core. Each Visium spot gets estimated proportions of your Level 2 cell types.

```python
# cell2location for VISIUM only
import cell2location

# Train reference signature model from your scRNA-seq
cell2location.models.RegressionModel.setup_anndata(
    adata_ref, batch_key='dataset_id', labels_key='cell_type_level2'
)
mod_ref = cell2location.models.RegressionModel(adata_ref)
mod_ref.train(max_epochs=250, use_gpu=True)

# Apply to Visium spatial data
cell2location.models.Cell2location.setup_anndata(adata_visium, batch_key='sample')
mod_spatial = cell2location.models.Cell2location(
    adata_visium,
    cell_state_df=mod_ref.export_posterior(adata_ref)
)
mod_spatial.train(max_epochs=30000, use_gpu=True)

# Result: per-spot estimated abundance of each cell type
# Visualize: spatial heatmaps showing SPP1+ Mac gradient, Basaloid gradient, etc.
```

### Track B: Xenium Data (Direct Label Transfer)
**Datasets**: Vannan 2025 (GEO per paper), Mayr 2024 (Zenodo 10015169, Xenium subset)
**Method**: scVI/scANVI label transfer OR scanpy.tl.ingest
**Purpose**: Map your scRNA-seq cell type annotations directly onto individually segmented Xenium cells at true single-cell resolution. No deconvolution needed.

```python
import scvi

# Option 1: scANVI label transfer (preferred — probabilistic)
scvi.model.SCANVI.setup_anndata(adata_ref, labels_key='cell_type_level2', 
                                  unlabeled_category='Unknown')
model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category='Unknown')
model.train(max_epochs=20)

# Transfer labels to Xenium cells
predictions = model.predict(adata_xenium)
adata_xenium.obs['predicted_celltype'] = predictions

# Option 2: Simpler ingest-based transfer
sc.tl.ingest(adata_xenium, adata_ref, obs='cell_type_level2')

# Result: each individual Xenium cell gets a cell type label
# Now you can do TRUE single-cell spatial analyses:
# - Direct co-localization of SPP1+ Mac and KRT17+ Basaloid cells
# - Cell-cell distance measurements (see Amendment 2)
# - Neighborhood composition analysis
```

### Impact on Figure 6:
- **6a–b**: Visium (cell2location) → broad tissue architecture and regional gradients
- **6c–d**: Xenium (label transfer) → single-cell resolution co-localization
- **6e–f**: Xenium-based distance and neighborhood analyses (see Amendment 2)

---

## AMENDMENT 2: Distance-Dependent Spatial Signaling Gradient
### Status: ✅ ADOPTED (genuinely clever analysis; strengthens causality claim)

**The problem in v2:** Squidpy co-occurrence analysis (v2 Section 3.7.2) shows that SPP1+ macrophages and aberrant basaloid cells are found in the same spatial neighborhoods. But a reviewer can argue: "Co-localization doesn't prove communication — they might just both accumulate in areas of high collagen density."

**The v3 fix:** Prove a **distance-dependent signaling gradient** using Xenium single-cell coordinates. Show that as distance from SPP1+ macrophage to KRT8+ cell decreases, the IPF-fate TF activity (and its target gene expression) increases exponentially.

```python
from scipy.spatial import cKDTree
import numpy as np
import pandas as pd

# Get spatial coordinates for relevant cell types from Xenium data
spp1_mac = adata_xenium[adata_xenium.obs['predicted_celltype'] == 'SPP1_Mac']
krt8_cells = adata_xenium[adata_xenium.obs['predicted_celltype'] == 'Transitional_KRT8']

# Build KD-tree from macrophage positions
mac_coords = np.column_stack([
    spp1_mac.obs['x_centroid'].values,
    spp1_mac.obs['y_centroid'].values
])
tree = cKDTree(mac_coords)

# For each KRT8+ cell, find distance to nearest SPP1+ macrophage
krt8_coords = np.column_stack([
    krt8_cells.obs['x_centroid'].values,
    krt8_cells.obs['y_centroid'].values
])
distances, _ = tree.query(krt8_coords, k=1)

# Add distances to the KRT8+ cell metadata
krt8_cells.obs['dist_to_nearest_spp1_mac'] = distances

# Extract the Master TF activity score (from Decoupler, transferred to spatial)
# and key target gene expression
# Plot: x-axis = distance (µm), y-axis = TF activity or target gene expression
# Expected: exponential decay — closer to macrophage = higher TF activity

# Formal statistical test: Spearman correlation + distance bins
from scipy.stats import spearmanr
rho, pval = spearmanr(distances, krt8_cells.obs['SOX9_activity'])
print(f"Spearman rho = {rho:.3f}, p = {pval:.2e}")

# Bin by distance (0-25µm, 25-50µm, 50-100µm, 100-200µm, >200µm)
# and show mean TF activity per bin (barplot with error bars)
```

**Why this is powerful:** This creates a dose-response curve in physical space. If SOX9 (or whatever the master TF turns out to be) activity is 3× higher in KRT8+ cells within 25 µm of an SPP1+ macrophage vs. those >200 µm away, that's strong evidence of a spatially-mediated signaling effect. Combine with LIANA+ predictions of the specific ligand-receptor pair, and you have a multi-evidence causality chain.

### Impact on Figures:
- New panel in **Figure 6e**: Distance-decay scatter plot (distance to SPP1+ Mac vs. TF activity in KRT8+ cells)
- New panel in **Figure 6f**: Binned distance barplot showing mean TF activity at different distance ranges

---

## AMENDMENT 3: Epigenetic Chromatin Validation
### Status: ⚠️ ADOPTED WITH CORRECTED ACCESSION (Gemini hallucinated the GEO ID)

**The concept is correct:** Anchoring the Decoupler-predicted TF switch in physical chromatin accessibility data is a powerful upgrade. If the master TF's binding motif is in open chromatin at the IPF bifurcation but closed in healthy cells, that upgrades the finding from "computational inference" to "epigenetically validated."

**⚠️ CRITICAL FACT-CHECK: GSE286182 does not exist.**
Gemini cited "Kaminski Lab 2025 snATAC-seq dataset (GSE286182)" — this accession number returns no results on GEO. It appears to be hallucinated.

**The verified IPF snATAC-seq dataset is:**
- **Valenzi et al. 2023** — GSE214085
- 3 IPF + 2 control lungs, snATAC-seq + scRNA-seq
- Published in *European Respiratory Journal*
- Confirmed open-access on GEO
- Already identified in the v1 dataset inventory

**Additionally**, the LungMAP consortium (Wang et al. 2021, eLife) has **snATAC-seq from healthy human lungs** (part of the LungMAP datasets), providing the baseline chromatin landscape.

### The v3 Implementation (Session 5b):

```python
# Download Valenzi et al. 2023 snATAC-seq data
# GEO: GSE214085

import muon as mu
import scanpy as sc

# Load snATAC-seq data
adata_atac = sc.read_h5ad('valenzi_2023_snATAC.h5ad')

# Focus on fibroblast and epithelial compartments
# (Valenzi focused on fibroblasts — epithelial coverage may be limited)

# Key analysis: Motif accessibility for the Master TF
# 1. Compute TF motif deviation scores (chromVAR)
# 2. Compare motif accessibility: IPF vs. Control
# 3. Specifically check: Is the SOX9/TP63/HIF1A motif 
#    differentially accessible in IPF?

# Using ArchR (via rpy2) or Signac (via rpy2) or pure Python with episcanpy
import episcanpy as epi

# Compute motif enrichment
epi.tl.motif_enrichment(
    adata_atac,
    motif_db='JASPAR2022',
    groupby='disease'  # IPF vs. Control
)

# Extract chromVAR deviation scores for candidate TFs
# Plot: Chromatin accessibility at SOX9 binding sites in IPF vs. Control
```

### Important caveats to document in Methods:
1. The Valenzi dataset has only 3 IPF + 2 control samples (small n for ATAC)
2. It focused primarily on fibroblasts (TWIST1 was the key finding) — epithelial cell coverage for KRT8+ cells may be limited
3. Frame as "supporting evidence from an independent epigenetic modality" not as definitive proof
4. If epithelial ATAC data is insufficient, the fibroblast TWIST1 finding still validates the broader concept that disease-specific TFs have corresponding chromatin accessibility changes

### Impact on Figures:
- New **Extended Data Figure**: Chromatin accessibility at Master TF motif sites in IPF vs. Control (Valenzi data)
- Reference in **Figure 3** panel or legend: "Consistent with open chromatin at [TF] binding sites in IPF lung (Extended Data Fig. X)"

---

## AMENDMENT 4: Non-COVID ARDS Bias Defense
### Status: ⚠️ PARTIALLY ADOPTED (concept is right, but execution requires caution)

**The concept is correct:** Using COVID-19 as the sole proxy for "acute lung injury" is a vulnerability. A clinician reviewer will ask whether the ARDS bifurcation is universal or SARS-CoV-2–specific.

**⚠️ FACT-CHECK on Gemini's suggestion:**
Gemini suggested GSE151263 (Jiang et al. 2020) as a "non-COVID ARDS" dataset. This dataset exists and IS publicly available, BUT:
- It contains **PBMCs** (peripheral blood), **NOT lung tissue**
- It profiles 3 sepsis+ARDS and 4 sepsis-only patients
- There are **NO epithelial cells** — you cannot map KRT8+ transitional cells from blood
- It is useful for monocyte/macrophage comparisons but NOT for the epithelial bifurcation analysis

**The HLCA Census likely does NOT contain labeled "bacterial pneumonia" or "sepsis" cells.** The HLCA focuses on lung tissue samples; non-COVID ARDS lung tissue scRNA-seq data effectively does not exist in public repositories (confirmed in our v1 dataset inventory).

### The Realistic v3 Strategy (3-tiered defense):

**Tier 1: Leverage Grant et al. 2021 (GSE155249)**
This is the NU SCRIPT dataset that profiled BALF from 10 severe COVID patients AND included non-COVID pneumonia comparisons (bacterial/viral). While the non-COVID arm used flow cytometry and bulk RNA-seq rather than scRNA-seq, you can still use the available scRNA-seq from both COVID and non-COVID BALF:

```python
# Download Grant et al. 2021 BALF scRNA-seq
# GSE155249 — includes non-COVID pneumonia patients

# Check if non-COVID BALF cells are present in scRNA-seq portion
# Grant profiled 10 COVID patients by scRNA-seq
# The non-COVID comparison was primarily flow/bulk
# BUT: Check if any non-COVID scRNA-seq cells exist in the deposit
```

**Tier 2: Check CellxGene Census for non-COVID annotations**
During Session 1, explicitly query:
```python
# Check what disease labels exist in the HLCA census
census = cellxgene_census.open_soma()
# Query all unique disease values
disease_values = (
    census["census_data"]["homo_sapiens"]
    .obs.read(column_names=["disease"])
    .concat()
    .to_pandas()["disease"]
    .unique()
)
print("Available disease labels:", sorted(disease_values))
# Look for: 'pneumonia', 'bacterial pneumonia', 'sepsis', 
# 'influenza', 'acute respiratory distress syndrome'
```

**Tier 3: Transparent Discussion Section (the honest approach)**
If non-COVID ARDS single-cell data simply isn't available (which is likely), write this in the Discussion:

> "A limitation of this study is that our acute lung injury arm is derived entirely from COVID-19 patients. While the KRT8+ transitional state and SPP1+ macrophage programs have been independently described in non-COVID contexts [cite mouse BLM models, pre-COVID studies], direct single-cell validation in human non-COVID ARDS lung tissue awaits future studies. We note that Grant et al. (2021) showed that macrophage–T cell feedback circuits in COVID-ARDS parallel those in bacterial pneumonia, supporting generalizability. Our findings generate a testable prediction: that non-COVID ARDS KRT8+ cells should show the same inflammatory/apoptotic branch rather than the basaloid branch observed in IPF."

This framing turns the limitation into a future direction — and a testable prediction is actually a strength.

### Impact: No figure changes. Adds to Session 1 (Census metadata check) and Discussion paragraph 5.

---

## AMENDMENT 5: Redefine "Apoptotic" Terminal State Semantics
### Status: ✅ ADOPTED (this is an important technical precision issue)

**The problem:** scRNA-seq QC routinely filters cells with high mitochondrial reads and low UMI counts — which are the hallmarks of dying/apoptotic cells. If you label the ARDS terminal state "Apoptotic," a reviewer can argue: "You filtered out the truly apoptotic cells in QC. What you're calling 'apoptotic' is just low-quality cells that escaped your filters."

**The v3 fix:** Rename the ARDS terminal state from "Apoptotic" to **"Inflammatory Arrest"** or **"Senescent-Inflammatory"** throughout the paper.

### The biological rationale:
The state you're actually capturing in scRNA-seq is NOT mid-apoptosis (those cells are already dead/fragmented). What survives QC is the **pre-apoptotic inflammatory arrest state** — cells that are:
- **Cell-cycle arrested**: CDKN1A (p21), CDKN2A (p16) — senescence markers
- **Actively producing inflammatory cytokines**: IL1B, CXCL8, IL6, CCL2
- **Responding to IFN signaling**: ISG15, IFITM3, MX1, OAS1
- **Expressing stress/damage markers**: GADD45A/B, ATF3, DDIT3
- **NOT yet executing apoptosis**: low CASP3 activation (which requires protein-level detection anyway — mRNA for caspases is unreliable)

### Updated terminology throughout the paper:

| v2 term | v3 term | Rationale |
|---------|---------|-----------|
| "Apoptotic terminal state" | **"Inflammatory arrest"** | What scRNA-seq actually captures |
| "Inflammatory apoptosis" | **"Senescent-inflammatory fate"** | More accurate at transcript level |
| "CASP3/8 markers" | **"CDKN1A/2A, IL1B, ISG15, GADD45A"** | Verifiable by scRNA-seq |
| "Cell death branch" | **"Inflammatory arrest branch"** | Scientifically defensible |

### Updated CellRank terminal states (for Session 4):
```python
# v2 (problematic):
# g.set_terminal_states(["AT1_mature", "Aberrant_Basaloid", "Apoptotic"])

# v3 (corrected):
g.set_terminal_states([
    "AT1_mature",           # Successful repair
    "Aberrant_Basaloid",    # IPF fibrotic fate (KRT17+, COL1A1+)
    "Inflammatory_Arrest"   # ARDS inflammatory fate (CDKN1A+, ISG15+, IL1B+)
])

# Define Inflammatory_Arrest by markers:
inflammatory_arrest_markers = [
    'CDKN1A', 'CDKN2A',    # Senescence / cell cycle arrest
    'IL1B', 'CXCL8', 'IL6', # Inflammatory cytokines
    'ISG15', 'IFITM3', 'MX1', # IFN response
    'GADD45A', 'ATF3',      # Stress response
    'CCL2', 'CXCL10'        # Chemokine production
]
```

### Impact: Terminology change throughout the paper. No structural changes to figures.

---

## ADDITIONAL EXECUTION TIPS ADOPTED FROM GEMINI

### A. Use CL Ontology IDs for Census Queries (Session 1)

Cell type string names can change between HLCA versions. Use Cell Ontology (CL) IDs for robust filtering:

```python
# Instead of: obs_value_filter = "cell_type == 'type II pneumocyte'"
# Use: obs_value_filter = "cell_type_ontology_term_id == 'CL:0002063'"

# Key CL IDs for our target cell types:
# AT2: CL:0002063
# AT1: CL:0002062  
# Alveolar macrophage: CL:0000583
# Classical monocyte: CL:0000860
# Fibroblast: CL:0000057
# Note: KRT8+ transitional, aberrant basaloid — these may NOT have CL IDs
# (they're disease-associated states, not standard cell types)
# → Fall back to string matching or marker-based identification
```

### B. Covariate Control in All Statistical Tests

IPF patients skew older, male, and often have smoking history. COVID patients are more demographically diverse. Every comparison must control for this:

```python
# MiloR design formula (Session 3):
milo.da_nhoods(mdata, design="~ disease + sex + age_group")

# Decoupler differential TF activity (Session 5):
# Use pseudobulk aggregation per donor, then DESeq2 with covariates
# formula: ~ disease + sex + age_group

# CellRank (Session 4):
# Document donor demographics per disease group in Table 1
# Show that bifurcation findings hold when stratified by age/sex
```

### C. Memory Management for Biomni Sessions

Add to EVERY Claude Code prompt for Biomni sessions:

```python
# Standard memory-safe imports and practices
import gc
from scipy.sparse import issparse, csr_matrix

# After every major operation:
del intermediate_variable
gc.collect()

# Keep matrices sparse
if not issparse(adata.X):
    adata.X = csr_matrix(adata.X)

# For subsetting: always .copy() to prevent views holding full data
adata_sub = adata[mask].copy()
```

---

## UPDATED SESSION PLAN (additions only)

The v2 session plan (Sessions 1–12) remains unchanged except:

### Session 5b (NEW — Day 16-17): Epigenetic Validation
```
"Download the Valenzi et al. 2023 snATAC-seq data (GSE214085).
Focus on fibroblast and any available epithelial cells.
Using episcanpy or chromVAR (via rpy2):
1. Compute TF motif deviation scores
2. Compare motif accessibility: IPF vs. Control
3. Specifically test whether our master switch TF motif is 
   differentially accessible
4. Generate chromatin accessibility heatmap for Extended Data"
```

### Session 9 UPDATE: Split Spatial Pipeline
```
"Session 9a — Visium Track (cell2location):
Train reference model on scRNA-seq. Deconvolve Franzén IPF Visium 
and Mothes COVID Visium. Generate spatial abundance maps.

Session 9b — Xenium Track (label transfer):
Transfer cell type labels to Vannan Xenium data using scANVI.
Run distance-dependent gradient analysis: distance from KRT8+ cells 
to nearest SPP1+ Mac vs. TF activity.
Use scipy.spatial.cKDTree. Generate scatter plots and binned barplots."
```

---

## v3 SUMMARY TABLE

| Amendment | Source | Status | Impact |
|-----------|--------|--------|--------|
| 1. Split Visium/Xenium pipeline | Gemini | ✅ Adopted | Prevents methodological error; strengthens spatial analysis |
| 2. Distance-dependent gradient | Gemini | ✅ Adopted | Adds spatial causality evidence to Figure 6 |
| 3. snATAC chromatin validation | Gemini (concept) | ⚠️ Adopted, accession corrected | Upgrades TF finding; uses GSE214085 (not hallucinated GSE286182) |
| 4. Non-COVID ARDS defense | Gemini | ⚠️ Partially adopted | Census check + Discussion framing; Gemini's GSE151263 is PBMCs not lung |
| 5. "Apoptotic" → "Inflammatory Arrest" | Gemini | ✅ Adopted | Prevents QC artifact criticism |
| A. CL Ontology IDs | Gemini | ✅ Adopted | More robust Census queries |
| B. Covariate control | Gemini | ✅ Adopted | Prevents confounding criticism |
| C. Memory management | Gemini | ✅ Adopted | Practical Biomni stability |

---

## THE COMPLETE DOCUMENT SET FOR CLAUDE CODE SESSIONS

When you start working, reference these documents in order:
1. **v2 Research Plan** — The master structure (figures, narrative, sessions)
2. **v3 Amendments** (this document) — Patches to apply on top of v2
3. **v1 Dataset Inventory** — Verified accessions and links for all datasets

You're ready for Session 1. The first deliverable remains: **verify that KRT8+ transitional cells exist in both COVID and IPF within the HLCA Census.**
