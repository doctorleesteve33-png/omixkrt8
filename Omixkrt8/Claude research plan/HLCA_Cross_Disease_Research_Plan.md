# Research Plan: Cross-Disease Single-Cell Atlas Analysis of Lung Injury
## Convergent Cellular Programs of Aberrant Repair in COVID-ARDS and IPF

**Target Journal Tier**: Nature Medicine (primary) → Cell (backup) → Nature Communications (safety)
**Estimated Timeline**: 4–6 months (computational + writing)
**Tools**: CZ CELLxGENE HLCA, Claude Code, Biomni Lab

---

## 1. SCIENTIFIC RATIONALE & HYPOTHESIS

### 1.1 The Knowledge Gap

COVID-ARDS and IPF represent the archetypal acute vs. chronic lung injury paradigms, yet emerging evidence hints at deep mechanistic convergence:

- Both diseases feature **SPP1+ profibrotic macrophages** (Sikkema et al. 2023 HLCA; Morse et al. 2019)
- Both produce **KRT5−/KRT17+ aberrant basaloid epithelial cells** (Adams et al. 2020; Melms et al. 2021)
- Both activate **CTHRC1+ pathologic fibroblasts** (Tsukui et al. 2020; Bharat et al. 2020)
- Post-COVID pulmonary fibrosis clinically resembles IPF (radiologically and functionally)

**However, no study has systematically compared these programs at single-cell resolution using harmonized, batch-corrected data.** Previous comparisons relied on visual inspection of separate UMAPs from independently processed datasets, making quantitative cross-disease inference unreliable.

### 1.2 Central Hypothesis

> Acute (COVID-ARDS) and chronic (IPF) lung injury converge on shared "aberrant repair modules" — coordinated gene programs operating across macrophage, epithelial, and fibroblast compartments — that represent conserved injury-response failure. These modules share master transcriptional regulators and intercellular communication circuits, but diverge in their upstream triggers and temporal dynamics, explaining why COVID-ARDS can resolve while IPF cannot.

### 1.3 Why This Paper Will Be High-Impact

| Factor | Justification |
|--------|--------------|
| **Clinical urgency** | Post-COVID fibrosis affects millions; understanding whether it follows IPF biology determines treatment strategy (pirfenidone/nintedanib vs. immunomodulation) |
| **Methodological rigor** | First study to leverage HLCA harmonization for cross-disease comparison, eliminating the batch-effect confounding that plagued earlier meta-analyses |
| **Scale** | >400K cells from two diseases in a unified embedding, far exceeding any prior cross-disease comparison |
| **Actionable output** | Shared regulators/targets → drug repurposing candidates; divergent features → precision diagnostics |
| **Validation** | Spatial transcriptomics validation using independent public Visium/Xenium data |

---

## 2. DATA ACQUISITION STRATEGY

### 2.1 Primary Data: CZ CELLxGENE HLCA

**Do NOT download raw .fastq or uncorrected .h5ad from GEO.** Instead:

```
Source: CZ CELLxGENE Discover Portal
URL: https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293
Format: .h5ad (AnnData) — harmonized, batch-corrected, cell-ontology-annotated
```

**Core datasets to extract from HLCA:**

| Dataset | GEO | Disease | Cells (approx.) | Why Include |
|---------|-----|---------|-----------------|-------------|
| Adams et al. 2020 | GSE136831 | IPF + controls | ~312,000 | Largest IPF scRNA-seq; aberrant basaloid cells |
| Habermann et al. 2020 | GSE135893 | IPF/PF + controls | ~114,000 | Independent IPF cohort; KRT17+ cells |
| Melms et al. 2021 | GSE171524 | COVID-ARDS (lethal) | ~116,000 | Benchmark COVID lung autopsy |
| Reyfman et al. 2019 | GSE122960 | IPF + controls | ~76,000 | Northwestern cohort; overlaps with Grant |
| Additional HLCA healthy | Various | Healthy controls | Variable | Establish normal baseline |

**Total: ~600,000+ cells, 3 disease groups (COVID-ARDS, IPF, Healthy), unified cell ontology.**

### 2.2 HLCA Download Protocol (Claude Code)

```python
# Step 1: Install cellxgene-census
# pip install cellxgene-census

import cellxgene_census
import scanpy as sc

# Step 2: Open the HLCA census and query specific datasets
# Option A: Direct .h5ad download from CellxGene collection page
# Option B: Use Census API for programmatic access

# Key fields in HLCA to leverage:
# - obs['cell_type_ontology_term_id']  → standardized CL ontology
# - obs['disease_ontology_term_id']    → MONDO terms
# - obs['tissue_ontology_term_id']     → UBERON terms
# - obs['suspension_type']             → cell vs. nucleus
# - obsm['X_scVI']                     → batch-corrected latent space (scVI)
# - raw.X                              → raw counts for DE analysis

# Step 3: Subset to relevant diseases
# COVID-19: MONDO:0100096
# Idiopathic pulmonary fibrosis: MONDO:0008345
# Normal/healthy: filter by 'normal' in disease field
```

### 2.3 Spatial Validation Data (Downloaded Separately)

These are NOT in HLCA and must be fetched from their original repositories:

| Study | Platform | Data | URL |
|-------|----------|------|-----|
| Franzén et al. 2024 | 10x Visium | IPF spatial | ArrayExpress E-MTAB-14121 |
| Mayr et al. 2024 | Visium + Xenium | IPF spatial niches | Zenodo 10012934 / 10015169 |
| Vannan et al. 2025 | Xenium + Visium HD | IPF spatial (largest) | GEO (per paper) |
| Mothes et al. 2023 | 10x Visium | COVID lung spatial | Zenodo 7533977 |
| Delorey et al. 2021 | GeoMx DSP | COVID multi-organ spatial | GEO GSE171668 |

---

## 3. COMPUTATIONAL ANALYSIS PIPELINE

### Phase 1: Data Preparation & Quality Control (Week 1–2)

**3.1.1 HLCA Subset Extraction**
- Extract COVID-ARDS, IPF, and matched healthy control cells from HLCA
- Retain all metadata fields (donor age, sex, smoking status, sample site)
- Keep BOTH the scVI latent representation (for clustering/embedding) AND raw counts (for DE/GRN)

**3.1.2 Cell Type Harmonization Audit**
- Map HLCA cell ontology terms to a simplified hierarchy for this study
- Create a custom 3-level annotation:
  - Level 1: Compartment (Epithelial / Immune / Stromal / Endothelial)
  - Level 2: Major cell type (AT2, AT1, Macrophage, Fibroblast, etc.)
  - Level 3: Subtype/state (SPP1+ Mac, FABP4+ Mac, Aberrant basaloid, etc.)
- Verify that aberrant basaloid cells, profibrotic macrophages, and pathologic fibroblasts are annotated consistently across COVID and IPF datasets

**3.1.3 Batch Effect Assessment**
- Quantify residual batch effects using kBET, LISI, and silhouette scores
- If needed, perform additional Harmony or scVI correction within disease groups
- Generate integrated UMAP colored by disease, dataset, donor, and cell type

### Phase 2: Cross-Disease Differential Abundance (Week 2–3)

**3.2.1 Compositional Analysis (Figures 1–2)**
- Use **MiloR** (Dann et al. 2022, Nature Biotechnology) for differential abundance testing in neighborhoods
  - Compare: COVID-ARDS vs. Healthy, IPF vs. Healthy, COVID-ARDS vs. IPF
  - Control for: donor age, sex, smoking status
- Use **scCODA** (Büttner et al. 2021) as a complementary Bayesian compositional analysis
- Key question: Which cell states are **shared** between COVID-ARDS and IPF (expanded in both vs. healthy) vs. **disease-specific**?

**3.2.2 Expected Key Findings (based on literature)**
- Shared expansions: SPP1+ macrophages, aberrant basaloid cells, CTHRC1+ fibroblasts
- COVID-specific: Inflammatory monocytes, IFN-stimulated macrophages, TP63+ basal progenitors
- IPF-specific: HAS1-high fibroblasts, senescent AT2 cells, proliferating SPP1/MERTK+ macrophages

### Phase 3: Consensus Meta-Program Discovery (Week 3–5) ★ CORE NOVELTY

**3.3.1 Consensus Non-negative Matrix Factorization (cNMF)**
This is the analytical centerpiece. Use **cNMF** (Kotliar et al. 2019, eLife) to discover gene expression programs (meta-programs) that operate across cell types and diseases WITHOUT prior assumptions.

```
For each major cell type (Macrophage, Epithelial, Fibroblast):
  1. Run cNMF on pooled COVID + IPF + Healthy cells
  2. Identify k programs (k=10-30, optimize via stability analysis)
  3. Score each cell for each program (usage matrix)
  4. Classify programs as:
     a. SHARED: active in both COVID-ARDS and IPF, low in healthy
     b. COVID-SPECIFIC: active in COVID only
     c. IPF-SPECIFIC: active in IPF only
     d. HOMEOSTATIC: active in healthy, reduced in disease
```

**3.3.2 Program Characterization**
- For each shared/specific program, extract:
  - Top 50 genes (by NMF loading weight)
  - GO/KEGG/Reactome enrichment
  - TF enrichment (via SCENIC regulons, see below)
  - Drug target overlap (via DGIdb, Open Targets)
- Name programs by their dominant functional signature (e.g., "Fibrotic ECM Remodeling Program", "IFN-Inflammatory Program", "Aberrant Differentiation Program")

**3.3.3 Cross-Compartment Program Correlation**
- Test whether activation of Program X in macrophages correlates with Program Y in fibroblasts within the same patient
- This reveals coordinated multi-compartment programs (the "aberrant repair module" hypothesis)
- Use patient-level program scores and Spearman correlation / mixed-effects models

### Phase 4: Gene Regulatory Network (GRN) Analysis (Week 4–6)

**3.4.1 pySCENIC (Aibar et al. 2017)**
- Run SCENIC on each major cell type separately (Macrophage, Epithelial, Fibroblast)
- Identify active regulons (TF + target gene sets) in each disease state
- Specifically interrogate:
  - **SPP1+ macrophage regulons**: PPARG, MAFB, CEBPB, NR1H3
  - **Aberrant basaloid regulons**: TP63, SOX9, KLF5, YAP1/WWTR1
  - **Pathologic fibroblast regulons**: TWIST1 (Valenzi 2023), FOXO3, RUNX2
- Identify **shared master regulators**: TFs whose regulons are active in both COVID-ARDS and IPF but not healthy

**3.4.2 CellOracle / GENIE3 (Complementary)**
- Use CellOracle (Kamimoto et al. 2023) for in silico TF perturbation
- Simulate: "What happens if you knock out TWIST1 in IPF fibroblasts? Does the same effect occur in COVID fibroblasts?"
- This provides causal directionality to the correlative SCENIC results

### Phase 5: Cell-Cell Communication Analysis (Week 5–7)

**3.5.1 CellChat v2 (Jin et al. 2021)**
- Run CellChat on COVID-ARDS, IPF, and Healthy separately
- Compare signaling pathway activities across conditions:
  - SPP1 signaling (macrophage → fibroblast)
  - TGFβ signaling (fibroblast → epithelial)
  - WNT signaling (epithelial ↔ fibroblast)
  - CXCL signaling (immune recruitment)
  - NOTCH signaling (epithelial differentiation)
- Identify **conserved circuits** (active in both diseases) vs. **rewired circuits** (disease-specific)

**3.5.2 NicheNet (Browaeys et al. 2020)**
- Use NicheNet to predict which ligands from macrophages/fibroblasts drive the aberrant basaloid state in epithelial cells
- Compare ligand rankings between COVID-ARDS and IPF contexts
- Identify shared upstream signals that could be therapeutically targeted

### Phase 6: Trajectory & Differentiation Analysis (Week 6–8)

**3.6.1 RNA Velocity (scVelo / CellRank)**
- Apply scVelo (Bergen et al. 2020) to the raw counts (spliced/unspliced)
- Note: This requires going back to the original BAM files or using the raw GEO data for spliced/unspliced counts (HLCA processed data may not have this)
- Alternative: Use CellRank's pseudotime approach (Lange et al. 2022) which works with the scVI latent space

**3.6.2 Key Trajectories to Map**
- AT2 → Aberrant basaloid → AT1 (or dead end?) in IPF vs. COVID
- Monocyte → SPP1+ profibrotic macrophage (shared trajectory?)
- Fibroblast → CTHRC1+ pathologic myofibroblast
- Compare: Are cells "stuck" at the same point in both diseases? Or do they traverse different paths to similar endpoints?

**3.6.3 Palantir / DPT Pseudotime**
- Compute diffusion pseudotime for epithelial cells in each condition
- Test: Is the "aberrant basaloid" state a terminal attractor in IPF but a transient state in COVID?
- This addresses the clinical question of why COVID fibrosis can sometimes resolve

### Phase 7: Spatial Validation (Week 7–9) ★ CRITICAL FOR HIGH IMPACT

**3.7.1 Spatial Deconvolution**
- Apply **cell2location** (Kleshchevnikov et al. 2022) or **RCTD** to deconvolve the scRNA-seq-defined cell states into spatial transcriptomics data
- Map the cNMF-defined meta-programs onto spatial coordinates
- Key question: Do the "shared aberrant repair modules" co-localize spatially?

**3.7.2 IPF Spatial Validation (Franzén + Mayr + Vannan data)**
- Map SPP1+ macrophage program, aberrant basaloid program, and pathologic fibroblast program to Visium spots
- Show that these programs co-localize in fibrotic foci (Visium/Xenium)
- Compare program spatial distribution across fibrotic severity gradient (using Franzén's mild→severe tissue blocks)

**3.7.3 COVID-ARDS Spatial Validation (Mothes + Delorey data)**
- Map the same programs to COVID lung Visium sections
- Show that shared programs occupy analogous spatial niches in DAD lesions
- Demonstrate that COVID-specific programs (IFN response) occupy distinct niches

**3.7.4 Niche Interaction Mapping**
- Use **COMMOT** (Cang et al. 2023) or **stLearn** for spatially-resolved cell-cell communication
- Validate CellChat predictions: Do the predicted signaling interactions occur between spatially adjacent cells?

### Phase 8: Therapeutic Target Prioritization (Week 8–10)

**3.8.1 Drug-Gene Interaction Database (DGIdb)**
- Query shared master regulators and key ligand-receptor pairs against DGIdb
- Focus on: druggable targets with existing compounds (repurposing potential)

**3.8.2 Open Targets Platform**
- Cross-reference shared regulators with genetic evidence from GWAS for IPF (Allen et al. 2017) and COVID-19 severity (Pairo-Castineira et al. 2021)
- Targets supported by BOTH transcriptomic and genetic evidence are highest priority

**3.8.3 Connectivity Map (CMap) Analysis**
- Query the shared "aberrant repair" gene signature against CMap/L1000
- Identify compounds that reverse the signature
- Prioritize compounds already in clinical trials for either disease

---

## 4. FIGURE PLAN (8 Main Figures + Extended Data)

### Figure 1: Cross-Disease Atlas Overview & Compositional Shifts
- **1a**: Study design schematic (HLCA → extraction → analysis pipeline)
- **1b**: Integrated UMAP colored by disease (COVID-ARDS / IPF / Healthy), split view
- **1c**: Integrated UMAP colored by harmonized cell type (Level 2 annotation)
- **1d**: Cell type composition barplots by disease (stacked, per-patient)
- **1e**: MiloR differential abundance beeswarm plot: COVID vs. Healthy
- **1f**: MiloR differential abundance beeswarm plot: IPF vs. Healthy
- **1g**: Venn diagram / UpSet plot: shared vs. disease-specific expanded populations

### Figure 2: Shared "Aberrant Repair" Cell States
- **2a**: Focused UMAP of macrophage compartment, colored by subtype
- **2b**: SPP1+ macrophage marker expression (SPP1, TREM2, CD9, GPNMB) — violin/dotplot by disease
- **2c**: Focused UMAP of epithelial compartment
- **2d**: Aberrant basaloid marker expression (KRT17, KRT5−, CDKN1A, TP63, CLDN4, COL1A1) by disease
- **2e**: Focused UMAP of fibroblast compartment
- **2f**: CTHRC1+ pathologic fibroblast markers by disease
- **2g**: Quantification: fraction of each "aberrant" state per patient, compared across diseases (box plots with significance)

### Figure 3: Meta-Program Discovery via cNMF ★ Central Figure
- **3a**: cNMF program stability plot (K selection)
- **3b**: Program usage heatmap (programs × cells, grouped by disease and cell type)
- **3c**: Classification of programs: Shared / COVID-specific / IPF-specific / Homeostatic (colored matrix)
- **3d**: Top genes for 3–4 key shared programs (lollipop/barplots)
- **3e**: GO enrichment for shared programs
- **3f**: Patient-level program scores: cross-compartment correlation matrix (showing coordinated multi-cell-type modules)

### Figure 4: Shared and Divergent Gene Regulatory Networks
- **4a**: SCENIC regulon activity heatmap (selected TFs × disease states), macrophages
- **4b**: Same for epithelial cells
- **4c**: Same for fibroblasts
- **4d**: Shared master regulators network visualization (Cytoscape-style)
- **4e**: CellOracle in silico KO: effect of shared TF perturbation in COVID vs. IPF fibroblasts
- **4f**: Regulon specificity score (RSS) plots identifying disease-specific TFs

### Figure 5: Cell-Cell Communication Circuits
- **5a**: CellChat aggregated signaling flow — Healthy
- **5b**: CellChat aggregated signaling flow — COVID-ARDS
- **5c**: CellChat aggregated signaling flow — IPF
- **5d**: Information flow comparison: pathways ranked by activity in each condition
- **5e**: Specific circuit deep-dives: SPP1→integrin (Mac→Fib), TGFβ (Fib→Epi), WNT (Epi↔Fib)
- **5f**: NicheNet top ligand predictions for aberrant basaloid state induction — COVID vs. IPF

### Figure 6: Divergent Trajectories Toward Shared Endpoints
- **6a**: CellRank-based fate maps for epithelial cells (AT2 → basaloid/AT1)
- **6b**: Pseudotime distribution of aberrant basaloid cells in COVID vs. IPF
- **6c**: Gene expression dynamics along pseudotime: key marker genes
- **6d**: Monocyte → SPP1+ macrophage trajectory comparison
- **6e**: Terminal state probability analysis: is basaloid a terminal attractor in IPF but not COVID?

### Figure 7: Spatial Validation of Shared Programs ★ Validation Figure
- **7a**: IPF Visium (Franzén data): spatial mapping of shared "aberrant repair" program score
- **7b**: IPF Xenium (Vannan/Mayr data): single-cell resolution co-localization of SPP1+ Mac + basaloid + CTHRC1+ Fib in fibrotic foci
- **7c**: COVID Visium (Mothes data): spatial mapping of same shared program
- **7d**: Spatial co-localization quantification: Moran's I or co-occurrence statistics
- **7e**: Niche-specific communication: COMMOT-predicted signaling in fibrotic niche (spatial)
- **7f**: Comparison of spatial niche architecture between COVID DAD and IPF fibrotic focus

### Figure 8: Therapeutic Implications
- **8a**: Shared druggable targets summary (bubble plot: target × evidence source × druggability)
- **8b**: CMap analysis: top compounds reversing the shared aberrant repair signature
- **8c**: GWAS overlap: shared targets with genetic evidence from IPF and COVID severity GWAS
- **8d**: Proposed model: convergent vs. divergent pathways, therapeutic intervention points (schematic)

### Extended Data / Supplementary (≥10 figures)
- QC metrics, batch correction validation (kBET/LISI)
- Complete cNMF program catalog (all programs, all cell types)
- Full SCENIC regulon lists
- Full CellChat pathway comparisons
- Additional spatial validations
- Sensitivity analyses (e.g., removing one dataset, re-running core analyses)
- Donor demographic tables
- Code availability / reproducibility documentation

---

## 5. PAPER STRUCTURE & NARRATIVE ARC

### Title (options)
1. "Convergent cellular programs of aberrant repair across acute and chronic human lung injury"
2. "A unified landscape of lung repair failure in COVID-ARDS and IPF reveals shared therapeutic targets"
3. "Cross-disease single-cell atlas analysis identifies conserved aberrant repair modules in human lung injury"

### Abstract (~250 words)
**Background** → **Gap** → **Approach** → **Key findings (3-4)** → **Significance**

### Introduction (4 paragraphs)
1. COVID-ARDS and IPF represent acute vs. chronic lung injury paradigms with distinct clinical courses but emerging biological convergence
2. Single-cell genomics has revealed disease-associated cell states in each disease independently (SPP1+ Mac, basaloid, etc.)
3. **GAP**: No systematic comparison using harmonized data; batch effects have prevented quantitative cross-disease inference; unknown whether shared cell states reflect identical or convergent gene programs
4. Here we leverage the HLCA to perform the first unified cross-disease analysis of >600K cells, discovering shared "aberrant repair modules" and identifying therapeutic targets

### Results (8 sections, matching 8 figures)
1. Harmonized cross-disease atlas reveals shared and disease-specific compositional shifts
2. COVID-ARDS and IPF share three "aberrant" cell states with quantitatively similar transcriptomes
3. cNMF identifies conserved meta-programs operating across cell types and diseases
4. Shared master regulators drive convergent cell states through parallel regulatory circuits
5. Cell-cell communication analysis reveals conserved profibrotic signaling circuits
6. Trajectory analysis shows convergent endpoints via divergent differentiation paths
7. Spatial transcriptomics validates co-localization of shared programs in fibrotic niches
8. Shared aberrant repair modules nominate druggable therapeutic targets

### Discussion (6 paragraphs)
1. Summary: first harmonized cross-disease atlas analysis
2. Biological significance: "aberrant repair module" as a unifying concept
3. Clinical implications: why post-COVID fibrosis may respond to anti-fibrotics; biomarker potential
4. The "attractor" model: why the same endpoint is reversible in acute but not chronic injury
5. Limitations (see below)
6. Future directions

---

## 6. METHODS SPECIFICATION

### 6.1 Data Acquisition
- HLCA v1.0 from CZ CELLxGENE (cite Sikkema et al. 2023)
- Spatial data from ArrayExpress, Zenodo, GEO (cite each original study)
- All data publicly available; no ethics approval required for secondary analysis of de-identified public data

### 6.2 Key Software (all Python-based, runnable in Claude Code + Biomni)

| Tool | Version | Purpose |
|------|---------|---------|
| scanpy | ≥1.9 | Core scRNA-seq analysis |
| scvi-tools | ≥1.0 | Batch correction, latent space |
| MiloR (milopy) | latest | Differential abundance |
| scCODA | latest | Compositional analysis |
| cNMF | ≥1.4 | Meta-program discovery |
| pySCENIC | ≥0.12 | Gene regulatory networks |
| CellChat (via rpy2) | v2 | Cell-cell communication |
| CellRank | ≥2.0 | Fate mapping / trajectories |
| cell2location | ≥0.1 | Spatial deconvolution |
| squidpy | ≥1.3 | Spatial analysis |
| decoupleR | latest | TF activity inference |

### 6.3 Statistical Framework
- Differential expression: Wilcoxon rank-sum (scanpy) + pseudobulk DESeq2 (for donor-level replication)
- Compositional: scCODA (Bayesian) + MiloR (graph-based)
- Multiple testing: Benjamini-Hochberg FDR < 0.05
- Effect sizes: log2FC thresholds (|log2FC| > 0.5 for DE; > 0.25 for program scores)
- Reproducibility: all random seeds fixed; Docker container for environment

---

## 7. EXECUTION PLAN: CLAUDE CODE + BIOMNI LAB

### 7.1 Division of Labor

**Claude Code (local/API — best for):**
- Data wrangling and preprocessing scripts
- Custom plotting code (matplotlib/seaborn publication figures)
- Statistical testing scripts
- Documentation and code organization
- Rapid prototyping of analysis functions
- Writing and editing manuscript text

**Biomni Lab (cloud compute — best for):**
- GPU-accelerated scVI model training (batch correction)
- pySCENIC (requires significant RAM, ~64GB+ for 600K cells)
- cNMF optimization runs (parallelizable)
- cell2location spatial deconvolution (GPU-accelerated)
- CellRank trajectory computation
- Any analysis requiring >16GB RAM

### 7.2 Week-by-Week Execution Timeline

```
WEEK 1-2:   Data acquisition + QC + UMAP generation
WEEK 2-3:   Compositional analysis (MiloR, scCODA) → Figure 1
WEEK 3-5:   cNMF meta-program discovery → Figure 3 (central)
WEEK 4-6:   SCENIC + CellOracle GRN analysis → Figure 4
WEEK 5-7:   CellChat + NicheNet communication → Figure 5
WEEK 6-8:   CellRank trajectory analysis → Figure 6
WEEK 7-9:   Spatial deconvolution + validation → Figure 7
WEEK 8-10:  Drug target analysis + synthesis → Figure 8
WEEK 10-12: Figure polishing + manuscript writing
WEEK 12-14: Internal review + revision
WEEK 14-16: Submission preparation
```

### 7.3 Critical Path Items (Potential Bottlenecks)

| Risk | Mitigation |
|------|-----------|
| HLCA download/parsing issues | Test CellxGene Census API first; fallback to direct .h5ad download |
| scVI latent space doesn't separate diseases well | Re-train scVI with disease as a covariate; or use Harmony |
| cNMF instability at chosen K | Run over K=5–40 with 100 iterations each; use consensus clustering |
| SCENIC memory overflow | Subsample to 100K cells per run; use GRNBoost2 GPU mode on Biomni |
| Spliced/unspliced counts unavailable in HLCA | Skip RNA velocity; rely on CellRank with diffusion pseudotime instead |
| Spatial data format incompatibilities | Use SpatialData (scverse) as universal format |

---

## 8. ANTICIPATED LIMITATIONS & PRE-EMPTIVE STRATEGIES

### 8.1 Limitations to Acknowledge

1. **Sample source bias**: IPF data from end-stage explants; COVID from lethal autopsy → both represent severe/terminal disease, not early stages
   - *Mitigation*: Discuss this explicitly; note that shared programs in end-stage disease may differ from early disease

2. **No non-COVID ARDS comparison**: All acute injury data is COVID-specific
   - *Mitigation*: Frame paper as COVID-ARDS vs. IPF (not "all ARDS"); discuss in limitations; propose non-COVID ARDS analysis as future work

3. **Computational cross-disease comparison, no new experimental data**
   - *Mitigation*: This is common for atlas-scale papers (see Sikkema 2023, Suo 2022); spatial validation provides independent confirmation; propose experimental validation as follow-up

4. **Harmonization may mask real biological differences**
   - *Mitigation*: Show results are robust to different batch correction methods (Harmony vs. scVI vs. scANVI); show key findings replicate in per-dataset analyses

### 8.2 Reviewer Anticipation Strategy

| Likely reviewer concern | Pre-emptive response |
|------------------------|---------------------|
| "Batch effects confound disease comparisons" | Show kBET/LISI/silhouette; repeat key analyses with pseudobulk aggregation per donor |
| "How do you know shared states aren't just generic stress response?" | Show specificity: shared states are NOT general stress/apoptosis but specific fibrotic programs absent in healthy tissue; compare to other injury types if data available |
| "No experimental validation" | Spatial validation + genetic evidence (GWAS) + drug database convergence provide orthogonal support; frame as hypothesis-generating resource |
| "Why not include COPD/other ILDs?" | Address in discussion; COPD data is in HLCA and could be added in revision |
| "Sample size concerns for some comparisons" | Report donor-level (not cell-level) n for all tests; use mixed-effects models |

---

## 9. JOURNAL STRATEGY

### 9.1 Target: Nature Medicine
**Fit**: Cross-disease mechanism with therapeutic implications; computational atlas with clinical relevance; precedent: Sikkema et al. 2023 (HLCA itself was published here)
**Format**: Article; ~5,000 words main text; 8 main figures; 10+ extended data figures
**Turnaround**: ~4-8 weeks to first decision

### 9.2 Backup: Cell
**Fit**: Systems-level biological insight; Resource paper aspect
**Format**: Article; ~7,000 words; flexible figure count

### 9.3 Safety: Nature Communications
**Fit**: Solid cross-disease comparison; technically rigorous
**Format**: Article; ~5,000 words; 8 figures
**Advantage**: Higher acceptance rate; still IF >14

### 9.4 Key Differentiators vs. Existing Literature

| Existing work | What we add |
|---------------|-------------|
| Sikkema 2023 (HLCA) | They built the atlas; we USE it for cross-disease biology |
| Adams 2020 / Habermann 2020 | They characterized IPF individually; we compare to COVID quantitatively |
| Melms 2021 | They characterized COVID individually; we compare to IPF quantitatively |
| Bharat 2020 | Small-scale COVID-IPF comparison (5 patients); we do atlas-scale (>600K cells) |
| Various COVID reviews | Qualitative observations of similarity; we provide quantitative meta-programs |

---

## 10. DATA & CODE SHARING PLAN

All code will be deposited on GitHub with a DOI (Zenodo archive):
```
GitHub repository structure:
├── 01_data_acquisition/
│   ├── download_hlca.py
│   └── download_spatial.py
├── 02_preprocessing/
│   ├── subset_diseases.py
│   └── qc_metrics.py
├── 03_compositional/
│   ├── milor_analysis.py
│   └── sccoda_analysis.py
├── 04_cnmf/
│   ├── run_cnmf.py
│   └── classify_programs.py
├── 05_grn/
│   ├── scenic_pipeline.py
│   └── celloracle_perturbation.py
├── 06_communication/
│   ├── cellchat_analysis.R
│   └── nichenet_analysis.R
├── 07_trajectory/
│   └── cellrank_analysis.py
├── 08_spatial/
│   ├── cell2location_deconv.py
│   └── spatial_validation.py
├── 09_drug_targets/
│   └── target_prioritization.py
├── figures/
│   └── [publication-ready figure scripts]
├── environment.yml
└── README.md
```

Processed AnnData objects will be deposited on Zenodo or CellxGene for community re-use.

---

## 11. FIRST CONCRETE STEPS (This Week)

### Step 1: Access HLCA from CellxGene
```bash
# In Claude Code or Biomni terminal:
pip install cellxgene-census scanpy anndata
```
- Download the HLCA core + extended (disease) datasets
- Inspect metadata fields, confirm COVID and IPF annotations are present
- Estimate total cell count and memory requirements

### Step 2: Generate Proof-of-Concept UMAP
- Create a single UMAP showing COVID-ARDS vs. IPF vs. Healthy cells
- Color by cell type; generate split panels by disease
- This is your "Figure 1b/c draft" — if this looks clean, the project is viable

### Step 3: Quick Win — SPP1+ Macrophage Comparison
- Subset macrophages from all three conditions
- Run basic DE analysis: SPP1+ macrophages in COVID vs. IPF
- If they share >70% of top DEGs (vs. healthy), the convergence hypothesis is strongly supported
- If they diverge substantially, pivot to a "convergent endpoints, divergent programs" framing

### Step 4: Assess Computational Requirements
- Profile memory/compute needs for the full dataset
- Determine which analyses can run locally (Claude Code) vs. need Biomni
- Set up the GitHub repository structure

---

## SUMMARY: WHAT MAKES THIS A TOP-TIER PAPER

1. **Right question, right time**: Post-COVID fibrosis is a clinical reality; the field needs mechanistic framework
2. **Methodological innovation**: First use of HLCA harmonization for cross-disease mechanism discovery (not just atlas-building)
3. **Scale**: >600K cells, 2 diseases, spatial validation — exceeds all prior comparisons
4. **Biological depth**: Not just "what's different" but "what gene programs are shared, what regulates them, where are they spatially, and can we drug them"
5. **Actionable output**: Named drug targets + CMap compounds → directly translatable
6. **Narrative clarity**: "Convergent aberrant repair" is a memorable concept that unifies disparate observations
7. **Reproducibility**: All public data, all code shared, harmonized foundation
