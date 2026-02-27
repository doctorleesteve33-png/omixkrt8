# Research Plan v2: Cross-Disease Single-Cell Atlas of Lung Repair Failure
## "Convergent Injury, Divergent Fate: The Macrophage–Epithelial–Fibroblast Axis in Acute vs. Chronic Human Lung Failure"

**Target**: Nature Medicine (primary) → Sci. Transl. Med. (backup) → Nature Communications (safety)
**Timeline**: 4–5 months (computational + writing)
**Tools**: CZ CELLxGENE HLCA Census API, Claude Code, Biomni Lab (cloud GPU/RAM)

---

## WHAT CHANGED IN v2 (AND WHY)

This plan integrates strategic upgrades from a parallel Gemini deep-think review. Below is a transparent audit of what was adopted, what was rejected, and what was already covered. This matters because your Claude Code sessions will reference this document directly.

### ✅ Adopted from Gemini (genuine improvements)

| Change | Rationale |
|--------|-----------|
| **Sharper hypothesis: KRT8+ transitional cell bifurcation** | The v1 hypothesis ("shared aberrant repair modules") is too broad for a single paper. Focusing on *why KRT8+ cells die in ARDS but become KRT17+ basaloid in IPF* is a testable, mechanistic question with a clear narrative arc. The convergent modules become the *substrate*, the divergent bifurcation becomes the *story*. |
| **Census API streaming instead of bulk .h5ad download** | Critical for Biomni memory management. The full HLCA .h5ad is ~30–50 GB. The Census API lets you stream only Epithelial + Myeloid + Fibroblast cells from the three disease conditions, reducing to ~2–5 GB in RAM. Prevents OOM crashes. |
| **Decoupler as primary TF activity tool (SCENIC demoted to secondary)** | Decoupler uses pre-computed TF–target databases (DoRothEA, CollecTRI) and runs in minutes on 500K cells. SCENIC requires ~64 GB RAM and hours of GRNBoost2 compute. Use Decoupler for main figures; reserve SCENIC for extended data validation on subsets. |
| **LIANA+ as meta-CCC framework (replaces standalone CellChat)** | LIANA+ wraps CellPhoneDB, CellChat, NATMI, NicheNet, and others into a consensus framework. Methodologically stronger than any single tool. Consensus scoring reduces false positives — a known reviewer concern. |
| **Bulk RNA-seq survival validation (GSE47460 LGRC cohort)** | My v1 plan lacked clinical outcome validation. Adding Kaplan-Meier survival analysis using the IPF signature from your TF bifurcation analysis closes the "translational gap" that reviewers will demand. This is the difference between Nature Communications and Nature Medicine. |
| **Session-by-session Claude Code execution strategy** | Breaking the analysis into discrete Claude Code sessions (one per analysis phase) prevents context window overflow and keeps each session focused. |

### ❌ Rejected from Gemini (where I disagree)

| Suggestion | Why rejected |
|------------|-------------|
| **Exclude fibroblasts entirely, focus only on Mac–Epi axis** | This is biologically wrong. The macrophage–epithelial–fibroblast triad is the core biology of fibrotic repair. SPP1+ macrophages signal to BOTH epithelial and fibroblast compartments. CTHRC1+ fibroblasts are the effectors of ECM deposition. Excluding them leaves the story incomplete and gives reviewers a free attack vector: "You ignored the effector cells." We include fibroblasts as a third compartment. |
| **Only 5 main figures** | Nature Medicine typically accepts 6–8 main figures for articles of this scope. 5 figures compresses too much; the spatial validation alone needs a full figure, and the survival/drug analysis needs another. We use **7 main figures** — tight but comprehensive. |
| **"Avoid standard clustering entirely"** | Extreme. You still need UMAP + Leiden for context and orientation (Figure 1). The point is that standard clustering is *necessary but not sufficient*. Advanced analyses (CellRank, cNMF, Decoupler) provide the mechanistic depth on top of the foundation. |
| **Implication that scANVI embedding is pre-computed in HLCA** | Needs verification. The HLCA core atlas uses scANVI, but the disease-extended atlas may use scVI. We verify this programmatically in Session 1 before building on assumptions. |

### 🔄 Already in v1 (Gemini's suggestions that were already covered)

cNMF meta-program discovery (Gemini missed this entirely — it's actually the most novel analytical element), CellRank trajectory analysis, spatial deconvolution with cell2location, CMap drug repurposing, reviewer anticipation strategy, detailed GitHub repository structure, and week-by-week timeline.

---

## 1. THE REFINED SCIENTIFIC HYPOTHESIS

### 1.1 The Two-Layer Question

**Layer 1 — Convergence (the substrate):**
In both COVID-ARDS and IPF, AT1 cells are destroyed. AT2 progenitors attempt regeneration but become trapped in a KRT8+ transitional state. SPP1+ macrophages accumulate at injury sites. CTHRC1+ fibroblasts deposit excess collagen. These are the *shared aberrant repair modules*.

**Layer 2 — Divergence (the story):**
Despite starting from the same KRT8+ intermediate, these cells take opposite paths:
- **In ARDS**: KRT8+ cells undergo **inflammatory apoptosis** (CASP3/8, IL1B exposure, IFNγ signaling) → cell death → potential resolution if surviving AT2 can finish repair
- **In IPF**: KRT8+ cells undergo **aberrant basaloid differentiation** (KRT17+, TP63low, CDKN1A+, COL1A1+) → permanent scar-forming state → no resolution possible

**The testable question:** What specific transcription factors are activated (or silenced) at the KRT8+ bifurcation point, and which macrophage-derived signals force the cell toward each fate?

### 1.2 Why This Framing Beats v1

| Aspect | v1 ("Convergent modules") | v2 ("Convergent injury, divergent fate") |
|--------|---------------------------|------------------------------------------|
| Narrative | Descriptive atlas | Mechanistic fate decision |
| Central figure | cNMF heatmap (interesting but abstract) | CellRank bifurcation map with TF switches (dramatic and intuitive) |
| Clinical hook | "Shared targets exist" | "We found the switch that determines if your lung heals or scars" |
| Drug target | Shared signatures → CMap | Bifurcation TFs → specific druggable targets + survival validation |
| Reviewer appeal | Solid but incremental | "This explains post-COVID fibrosis at the cellular decision level" |

### 1.3 The "One Sentence for Nature Medicine"

> We identify the transcription factor switch at the KRT8+ epithelial bifurcation point that determines whether injured human lung resolves (ARDS) or scars permanently (IPF), nominate druggable master regulators, validate their spatial co-localization with profibrotic macrophages, and demonstrate their prognostic significance for patient survival.

---

## 2. DATA ACQUISITION (Census API Streaming)

### 2.1 Primary Data: HLCA via CellxGene Census API

**DO NOT download the full .h5ad file.** Stream only what you need:

```python
# SESSION 1: Data Extraction (Claude Code or Biomni Jupyter)
# ============================================================

import cellxgene_census
import scanpy as sc
import numpy as np

# Open the Census (latest HLCA release)
census = cellxgene_census.open_soma(census_version="2023-12-15")

# Define filters BEFORE downloading — this is the key memory trick
DISEASE_FILTER = (
    "disease == 'normal' or "
    "disease == 'COVID-19' or "
    "disease == 'idiopathic pulmonary fibrosis'"
)

# Cell type filter: Epithelial + Myeloid + Fibroblast lineages
# Using cell_type_ontology_term_id for precise CL ontology filtering
# OR use broader compartment labels if available in HLCA metadata
CELL_TYPE_FILTER = (
    "cell_type in ['type II pneumocyte', 'alveolar macrophage', "
    "'classical monocyte', 'non-classical monocyte', "
    "'lung macrophage', 'fibroblast', 'myofibroblast', "
    "'adventitial cell', 'pericyte', "
    "'basal cell', 'club cell', 'AT2 transitional', "
    "'type I pneumocyte']"
)
# NOTE: Exact field names/values depend on HLCA metadata schema.
# Session 1 Task 1: Inspect available obs columns and values first.

# Stream and convert to AnnData
adata = cellxgene_census.get_anndata(
    census,
    organism="Homo sapiens",
    obs_value_filter=DISEASE_FILTER,
    # Column selection to reduce memory
    obs_column_names=[
        "cell_type", "disease", "tissue",
        "donor_id", "sex", "development_stage",
        "suspension_type", "assay",
        "dataset_id", "cell_type_ontology_term_id"
    ],
)

# IMMEDIATELY save to disk to free memory
adata.write_h5ad("hlca_subset_raw.h5ad")
print(f"Downloaded: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"Disease distribution:\n{adata.obs['disease'].value_counts()}")
print(f"Cell type distribution:\n{adata.obs['cell_type'].value_counts()}")

census.close()
```

**Fallback if Census API is unavailable or slow:**
Download specific collection .h5ad files from the CellxGene web UI, then subset in Python. The Adams (GSE136831) and Melms (GSE171524) datasets have dedicated CellxGene pages.

### 2.2 Verify Pre-Computed Embeddings

```python
# SESSION 1 continued: Check what's already computed
print("Available obsm keys:", list(adata.obsm.keys()))
# Expected: 'X_scVI' or 'X_scANVI' or 'X_umap'
# If scVI/scANVI latent space exists → use it directly
# If not → we'll need to run scVI ourselves (Session 2, on Biomni GPU)

print("Available layers:", list(adata.layers.keys()))
# Need: raw counts in adata.raw.X or adata.layers['counts']
# Required for: DE analysis, SCENIC, cNMF, RNA velocity
```

### 2.3 Secondary Data (Downloaded Separately)

| Purpose | Dataset | Source | Format |
|---------|---------|--------|--------|
| **IPF spatial (Visium)** | Franzén et al. 2024 | ArrayExpress E-MTAB-14121 | Visium .h5 + images |
| **IPF spatial (Xenium)** | Vannan et al. 2025 | GEO (per paper) | Xenium .h5ad |
| **IPF spatial (Xenium+Visium)** | Mayr et al. 2024 | Zenodo 10012934 | .h5ad + counts |
| **COVID spatial (Visium)** | Mothes et al. 2023 | Zenodo 7533977 | Visium .h5ad |
| **IPF bulk survival** | LGRC cohort | GEO GSE47460 | Expression matrix + clinical |
| **IPF GWAS summary stats** | Allen et al. 2020 | Public | Summary statistics |

---

## 3. COMPUTATIONAL PIPELINE (7 Phases)

### Phase 1: QC, Embedding, and Atlas Orientation (Week 1–2)

**Session 1–2 in Claude Code / Biomni**

**3.1.1 Quality Control**
```
- Filter: mitochondrial % < 20, nGenes 200–8000, nCounts 500–50000
- Remove doublets (Scrublet, if not already flagged in HLCA)
- Verify donor-level sample sizes per disease group
- Log-normalize + HVG selection (if raw counts; skip if HLCA already processed)
```

**3.1.2 Embedding Strategy**
- **If HLCA provides scVI/scANVI latent space**: Use directly for UMAP, neighbors, clustering
- **If not**: Run scVI on Biomni GPU (latent_dim=30, n_epochs=400, batch_key='dataset_id')
- Generate UMAP with 3 split panels: Healthy | COVID-ARDS | IPF
- Leiden clustering at multiple resolutions (0.5, 1.0, 1.5)

**3.1.3 Cell Type Harmonization (3-tier scheme)**

```
Level 1 (Compartment): Epithelial | Myeloid | Stromal
Level 2 (Major type):
  Epithelial: AT2, AT1, Transitional (KRT8+), Aberrant Basaloid (KRT17+), 
              Basal, Club, Ciliated
  Myeloid:    Alveolar Mac (FABP4+), SPP1+ Mac, Inflammatory Mono, 
              cDC1, cDC2, pDC
  Stromal:    Fibroblast, CTHRC1+ Myofibroblast, HAS1+ Fib, 
              Pericyte, Smooth Muscle
Level 3 (Fine state): Data-driven subclusters within each Level 2 type
```

**Key verification at this stage:**
- Can you identify KRT8+ transitional cells in BOTH COVID-ARDS and IPF?
- Can you identify KRT17+/KRT5− aberrant basaloid cells primarily in IPF?
- Can you identify SPP1+/TREM2+ macrophages in both diseases?
- If YES → proceed. If NO → re-examine cell type annotations and thresholds.

### Phase 2: Compositional Analysis — What Expands Where (Week 2–3)

**Session 3 in Claude Code**

**3.2.1 MiloR / pertpy Differential Abundance**

```python
# Using pertpy (Python implementation of Milo)
import pertpy as pt

# Build Milo object from the scVI-based neighbor graph
milo = pt.tl.Milo()
mdata = milo.load(adata)

# Assign cells to overlapping neighborhoods
milo.make_nhoods(mdata, prop=0.1)

# Count cells per neighborhood per sample (donor-level)
milo.count_nhoods(mdata, sample_col="donor_id")

# Test differential abundance
# Contrast 1: COVID-ARDS vs. Healthy
# Contrast 2: IPF vs. Healthy  
# Contrast 3: COVID-ARDS vs. IPF (CRITICAL — the divergence)
milo.da_nhoods(mdata, design="~ disease + sex + age_group")
```

**3.2.2 Compositional Analysis with scCODA (Bayesian validation)**
- Run scCODA as an orthogonal check on MiloR results
- Report both methods in the paper — consensus strengthens claims

**3.2.3 Expected Core Finding for Figure 1:**
- **Shared expansions** (both diseases vs. healthy): SPP1+ Mac, KRT8+ transitional, CTHRC1+ myofibroblast
- **COVID-specific**: IFN-stimulated macrophages (ISG15+, MX1+), inflammatory monocytes
- **IPF-specific**: KRT17+ aberrant basaloid, HAS1-high fibroblasts, proliferating SPP1/MERTK+ macrophages
- **Shared depletions**: FABP4+ resident alveolar macrophages, mature AT1 cells

### Phase 3: Trajectory & Bifurcation Mapping — The Paper's Central Figure (Week 3–5) ★

**Session 4–5 on Biomni (compute-intensive)**

**3.3.1 CellRank 2 Fate Mapping**

This is the analytical heart of the paper. You're reconstructing the *decision point* where an injured AT2 cell chooses its fate.

```python
import cellrank as cr
import scvelo as scv

# Subset to epithelial cells only for trajectory analysis
adata_epi = adata[adata.obs['compartment'] == 'Epithelial'].copy()

# Option A: If spliced/unspliced counts available → RNA Velocity
# scv.pp.moments(adata_epi)
# scv.tl.velocity(adata_epi, mode='dynamical')
# vk = cr.kernels.VelocityKernel(adata_epi)

# Option B (more likely): Use CellRank 2's PseudotimeKernel
# with diffusion pseudotime rooted at healthy AT2 cells
cr.tl.terminal_states(adata_epi, cluster_key='cell_type_level2')

# Define root: healthy AT2 cells
# Define terminal states: 
#   1. Mature AT1 (successful repair)
#   2. KRT17+ Aberrant Basaloid (IPF fate)
#   3. Apoptotic/Inflammatory (ARDS fate)

pk = cr.kernels.PseudotimeKernel(adata_epi, time_key="dpt_pseudotime")
pk.compute_transition_matrix()

# Compute absorption probabilities toward each terminal state
g = cr.estimators.GPCCA(pk)
g.compute_macrostates(n_states=5)
g.set_terminal_states(["AT1_mature", "Aberrant_Basaloid", "Apoptotic"])
g.compute_fate_probabilities()

# KEY OUTPUT: For each KRT8+ transitional cell, what is its probability
# of becoming Basaloid (IPF) vs. Apoptotic (ARDS) vs. AT1 (healthy repair)?
```

**3.3.2 Identifying the Bifurcation Point**
- Use CellRank's `compute_lineage_drivers()` to find genes that differ between the Basaloid and Apoptotic fates
- Compute a "bifurcation score" per cell: `P(Basaloid) - P(Apoptotic)`
- The cells where this score ≈ 0 are AT the bifurcation point
- Cells from IPF patients should have bifurcation scores shifted toward Basaloid
- Cells from ARDS patients should shift toward Apoptotic

**3.3.3 Gene Dynamics Along the Bifurcation**
- Plot expression of key genes across pseudotime, split by fate:
  - Shared early activation: KRT8, CLDN4, CDKN1A (senescence)
  - ARDS-fate genes: CASP3, CASP8, FAS, TRAIL, ISG15, IFNG-response genes
  - IPF-fate genes: KRT17, TP63, SOX9, COL1A1, TGFBI, MMP7
- This produces the classic "branching gene expression" plots for Figure 2

### Phase 4: Transcription Factor Switches — The Mechanism (Week 4–6) ★

**Session 6 in Claude Code + Biomni**

**3.4.1 Decoupler TF Activity Inference (PRIMARY — fast, scalable)**

```python
import decoupler as dc

# Get TF-target database (DoRothEA regulons, confidence A-C)
net = dc.get_dorothea(organism='human', levels=['A', 'B', 'C'])

# Run multivariate linear model (MLM) for TF activity scoring
dc.run_mlm(
    mat=adata_epi,       # Can also run on mac and fib subsets
    net=net,
    source='source',     # TF name
    target='target',     # Target gene
    weight='weight',     # Interaction confidence
    verbose=True,
    use_raw=True
)

# Result: adata_epi.obsm['mlm_estimate'] → per-cell TF activity scores
# adata_epi.obsm['mlm_pvals'] → significance

# KEY ANALYSIS: Compare TF activities at the bifurcation point
# Subset to KRT8+ transitional cells only
krt8_cells = adata_epi[adata_epi.obs['cell_type_level2'] == 'Transitional_KRT8']

# Split by disease
krt8_ards = krt8_cells[krt8_cells.obs['disease'] == 'COVID-19']
krt8_ipf = krt8_cells[krt8_cells.obs['disease'] == 'idiopathic pulmonary fibrosis']

# Differential TF activity: IPF-KRT8+ vs. ARDS-KRT8+
# This reveals which TFs are differentially activated at the SAME 
# transitional state but in different disease contexts
```

**3.4.2 Candidate Master Switches (Hypothesis-Driven)**

| TF | Expected Role | Evidence |
|----|---------------|----------|
| **SOX9** | Locks cells in aberrant progenitor state; blocks AT1 maturation | Strunz et al. 2020; active in IPF basaloid |
| **TP63** | Basal cell program driver; aberrant in distal lung | Melms 2021; Adams 2020 |
| **HIF1A** | Hypoxia-driven fibrotic signaling | Fibrotic foci are hypoxic; Goodwin 2018 |
| **YAP1/WWTR1** | Hippo pathway; mechanosensing in stiff fibrotic ECM | Liu et al. 2015; Xu et al. 2016 |
| **TWIST1** | Myofibroblast activation | Valenzi 2023 (only open-access IPF snATAC-seq) |
| **NF-κB (RELA)** | Inflammatory cascade; expected in ARDS path | Pro-apoptotic in epithelial cells |
| **IRF1/STAT1** | IFN-γ response; expected in ARDS path | Drives inflammatory death |
| **CEBPB** | Macrophage polarization; SPP1+ program | Active in profibrotic macrophages |

**3.4.3 pySCENIC (SECONDARY — extended data validation)**
- Run on Biomni GPU on a 50K-cell subsample (Epithelial only)
- Validates Decoupler findings with an independent GRN inference method
- If SCENIC and Decoupler agree on the master switch → very strong evidence

**3.4.4 CellOracle In Silico Perturbation (mechanistic direction)**
- Simulate: "What happens to KRT8+ cell fate probability if we knock out SOX9?"
- If KO of the IPF-enriched TF shifts fate probability toward AT1 repair → causal evidence
- This is the "prediction" that spatial/clinical data will validate

### Phase 5: cNMF Meta-Program Discovery — The Unbiased Layer (Week 4–6)

**Session 7 on Biomni**

*This is retained from v1 because Gemini missed it entirely, and it's actually the most methodologically novel element.*

**3.5.1 Consensus NMF on Each Compartment**

```python
from cnmf import cNMF

# Run cNMF on each compartment separately
for compartment in ['Epithelial', 'Myeloid', 'Stromal']:
    adata_sub = adata[adata.obs['compartment'] == compartment]
    
    # Initialize cNMF
    cnmf_obj = cNMF(
        output_dir=f'cnmf_{compartment}',
        name=f'{compartment}_programs'
    )
    
    # Prepare (select overdispersed genes)
    cnmf_obj.prepare(
        counts_fn=f'{compartment}_counts.h5ad',
        components=np.arange(5, 41),  # Test K=5 to 40
        n_iter=100,
        seed=42
    )
    
    # Factorize (parallelizable on Biomni)
    cnmf_obj.factorize()
    
    # Combine and select optimal K
    cnmf_obj.combine()
    cnmf_obj.k_selection_plot()  # Visual stability assessment
```

**3.5.2 Program Classification**
For each discovered program, classify as:
- **SHARED**: Active in both COVID-ARDS and IPF, low in Healthy
- **COVID-SPECIFIC**: Active only in COVID-ARDS
- **IPF-SPECIFIC**: Active only in IPF
- **HOMEOSTATIC**: Active in Healthy, suppressed in disease

**3.5.3 Cross-Compartment Coordination**
- Compute patient-level program scores
- Test: When "Fibrotic ECM Program" is high in fibroblasts, is "SPP1+ Program" also high in macrophages from the same patient?
- This reveals the multi-cell-type coordination that defines the "aberrant repair module"

### Phase 6: Cell-Cell Communication — Who Commands the Bifurcation (Week 5–7) ★

**Session 8 in Claude Code + Biomni**

**3.6.1 LIANA+ Meta-Framework (PRIMARY)**

```python
import liana as li

# Run LIANA+ with multiple methods for consensus
li.mt.rank_aggregate(
    adata,
    groupby='cell_type_level2',
    resource_name='consensus',  # Combines CellPhoneDB + CellChat + others
    expr_prop=0.1,
    verbose=True,
    use_raw=True
)

# Extract results
liana_res = adata.uns['liana_res']

# KEY COMPARISONS:
# 1. SPP1+ Mac → KRT8+ Transitional: What ligands in each disease?
# 2. SPP1+ Mac → CTHRC1+ Fibroblast: Shared or divergent?
# 3. CTHRC1+ Fib → KRT8+ Transitional: ECM signals?
```

**3.6.2 Hypothesis-Driven Pathway Focus**

| Pathway | Expected in ARDS | Expected in IPF | Biological Rationale |
|---------|-----------------|-----------------|---------------------|
| **TGFB1→TGFBR** | Moderate | Very High (Mac→Epi, Mac→Fib) | Master fibrotic signal |
| **SPP1→Integrins** | Present | Dominant | SPP1+ Mac → fibroblast activation |
| **WNT→Frizzled** | Low | High | Aberrant differentiation |
| **IL1B/TNF→TNFR** | Very High | Low | Inflammatory death signal |
| **IFNG→IFNGR** | Very High | Low | Apoptotic/inflammatory cascade |
| **NOTCH→DLL/JAG** | Variable | Aberrant | Epithelial fate specification |
| **FGF→FGFR** | Recovery signal? | Absent? | Normal repair pathway |

**3.6.3 NicheNet Ligand Prioritization**
- For the KRT8+ bifurcation: Which macrophage ligands best predict the IPF-fate TF activation vs. ARDS-fate TF activation?
- This links Phase 4 (TFs) to Phase 6 (signals) mechanistically

### Phase 7: Spatial & Clinical Validation — The "Reviewer-Proof" Layer (Week 6–9) ★★

**Session 9–11 across Claude Code + Biomni**

**3.7.1 Spatial Deconvolution: Map Cell States onto Tissue (cell2location)**

```python
import cell2location
import squidpy as sq

# Reference: your scRNA-seq AnnData with Level 2 cell type annotations
# Query: Spatial transcriptomics AnnData (Visium or Xenium)

# For IPF: Use Franzén 2024 Visium (fibrotic gradient) + Vannan 2025 Xenium
# For COVID: Use Mothes 2023 Visium

# Train cell2location reference model (on Biomni GPU)
cell2location.models.RegressionModel.setup_anndata(
    adata_ref,
    batch_key='dataset_id',
    labels_key='cell_type_level2'
)
mod = cell2location.models.RegressionModel(adata_ref)
mod.train(max_epochs=250, use_gpu=True)

# Map to spatial data
cell2location.models.Cell2location.setup_anndata(
    adata_spatial,
    batch_key='sample'
)
mod_spatial = cell2location.models.Cell2location(
    adata_spatial,
    cell_state_df=mod.export_posterior(adata_ref)
)
mod_spatial.train(max_epochs=30000, use_gpu=True)
```

**3.7.2 Spatial Validation Questions (each becomes a panel in Figure 6)**

| Question | Analysis | Expected result |
|----------|----------|-----------------|
| Do SPP1+ Mac and KRT17+ basaloid co-localize? | Squidpy co-occurrence, Moran's I | Yes, within fibrotic foci (IPF) |
| Are KRT8+ transitional cells at the edge of fibrotic foci? | Spatial deconvolution mapping | Yes — they are the "frontier" |
| Is the IPF-fate TF (e.g., SOX9) spatially restricted to fibrotic regions? | Map Decoupler TF scores onto space | Active in foci, not normal parenchyma |
| Does the COVID spatial pattern differ? | Same analysis on Mothes data | SPP1+ Mac more diffuse; KRT8+ near DAD |
| Do LIANA+-predicted signaling pairs occur between adjacent cells? | COMMOT or stLearn | Yes — predicted interactions are spatially proximal |

**3.7.3 Clinical Survival Validation (Kaplan-Meier)**

```python
# NEW IN v2 — adopted from Gemini's suggestion
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
import GEOparse

# Download LGRC IPF cohort (GSE47460)
# This has bulk RNA-seq + clinical outcomes for ~100+ IPF patients

# Step 1: Extract the "Bifurcation TF Signature"
# Top 20-50 genes from your Decoupler analysis that distinguish
# IPF-fate from ARDS-fate at the KRT8+ bifurcation
ipf_fate_genes = [...]  # From Phase 4 results

# Step 2: Score each patient for the signature
# Use ssGSEA or mean z-score of signature genes
from gseapy import ssgsea
ssgsea_scores = ssgsea(
    data=bulk_expression_matrix,
    gene_sets={'IPF_Fate_Signature': ipf_fate_genes},
    outdir=None,
    no_plot=True
)

# Step 3: Split patients into high/low signature groups (median split)
patients['signature_high'] = patients['score'] > patients['score'].median()

# Step 4: Kaplan-Meier + Cox regression
kmf = KaplanMeierFitter()
# Plot survival curves for high vs. low signature
# Report: HR, 95% CI, log-rank p-value

# Step 5: Multivariate Cox model (adjust for age, sex, FVC%, DLCO%)
cph = CoxPHFitter()
cph.fit(patients[['signature_score', 'age', 'sex', 'FVC_pct', 'time', 'event']],
        duration_col='time', event_col='event')
```

**3.7.4 Drug Repurposing (CMap/LINCS L1000)**

```python
import gseapy as gp

# Query LINCS L1000 with your IPF-fate signature
results = gp.enrichr(
    gene_list=ipf_fate_genes,
    gene_sets='LINCS_L1000_Chem_Pert_Consensus_Sigs',
    outdir='lincs_results'
)

# Also query: DrugMatrix, DSigDB, LINCS_L1000_Ligand_Perturbations
# Look for compounds that REVERSE the IPF-fate signature
# Prioritize FDA-approved drugs → immediate repurposing potential
```

---

## 4. FIGURE PLAN (7 Main Figures)

### Figure 1: The Unified Cross-Disease Landscape
- **1a**: Study design schematic (HLCA Census → filter → analysis pipeline)
- **1b**: Integrated UMAP, split by disease (Healthy | COVID-ARDS | IPF)
- **1c**: UMAP colored by Level 2 cell type
- **1d**: MiloR beeswarm: differential abundance across all three comparisons
- **1e**: Stacked barplot: cell type proportions per patient, grouped by disease
- **1f**: Venn/UpSet: shared vs. disease-specific expanded populations

### Figure 2: The KRT8+ Bifurcation — Divergent Fates from a Common Origin ★ CENTRAL
- **2a**: CellRank fate map (epithelial cells): AT2 → KRT8+ → (AT1 | Basaloid | Apoptotic)
- **2b**: Bifurcation score distribution: COVID vs. IPF patients (violin/ridge plot)
- **2c**: Fate probability heatmap per patient (rows=patients, columns=terminal states)
- **2d**: Gene expression dynamics along the two branches (branching lineplot)
  - ARDS branch: CASP3, FAS, ISG15, IFITM3 rising
  - IPF branch: KRT17, COL1A1, MMP7, SOX9 rising
  - Shared early: KRT8, CLDN4, CDKN1A rising then splitting
- **2e**: CellRank lineage driver genes: volcano plot of ARDS-fate vs. IPF-fate drivers

### Figure 3: Transcription Factor Switches at the Decision Point ★ MECHANISM
- **3a**: Decoupler TF activity heatmap: Top 30 TFs × cell states (Healthy AT2 → KRT8+ → Basaloid vs. Apoptotic)
- **3b**: TF activity along pseudotime, split by fate (line plots for top 6 TFs)
  - IPF-fate TFs: SOX9, TP63, YAP1 → rising on IPF branch
  - ARDS-fate TFs: STAT1, IRF1, NFKB1 → rising on ARDS branch
- **3c**: CellOracle in silico KO: SOX9 knockout shifts fate probability toward AT1
- **3d**: Regulon target network for the #1 master switch (network graph)
- **3e**: SCENIC validation (extended data reference): concordance with Decoupler

### Figure 4: cNMF Meta-Programs — Shared and Divergent Transcriptional Modules
- **4a**: cNMF K selection (stability plot)
- **4b**: Program × cell heatmap, grouped by disease and compartment
- **4c**: Shared vs. specific program classification matrix
- **4d**: Top genes for 3–4 key shared programs (bar plots)
- **4e**: Cross-compartment coordination: patient-level program correlation matrix
- **4f**: GO/Reactome enrichment for key programs

### Figure 5: The Macrophage Commands — Divergent Signals from SPP1+ Cells
- **5a**: LIANA+ consensus ligand-receptor dot plot: SPP1+ Mac → KRT8+ Transitional Epi
  - Side-by-side: COVID-ARDS vs. IPF contexts
- **5b**: LIANA+ consensus: SPP1+ Mac → CTHRC1+ Fibroblast
- **5c**: Information flow comparison: TGFβ, SPP1, WNT, TNF, IFNγ pathways across diseases
- **5d**: NicheNet: top ligands predicting IPF-fate TF activation vs. ARDS-fate TF activation
- **5e**: Proposed signaling model: which macrophage signals drive which fate

### Figure 6: Spatial Validation — The Fibrotic Niche Architecture ★ VALIDATION
- **6a**: IPF Visium (Franzén): cell2location mapping of SPP1+ Mac + Basaloid + CTHRC1+ Fib
- **6b**: IPF Xenium (Vannan): subcellular resolution co-localization in fibrotic foci
- **6c**: COVID Visium (Mothes): same cell states mapped, showing different spatial organization
- **6d**: Spatial TF activity: SOX9/TP63 activity mapped onto IPF tissue coordinates
- **6e**: Squidpy co-occurrence statistics: quantification of Mac–Basaloid spatial proximity
- **6f**: COMMOT spatially-resolved signaling: TGFβ signaling intensity maps

### Figure 7: Clinical Translation — Survival, Drugs, and the Therapeutic Window ★ IMPACT
- **7a**: Kaplan-Meier curves: IPF-fate signature (high vs. low) in LGRC cohort (GSE47460)
- **7b**: Forest plot: multivariate Cox regression (signature + clinical covariates)
- **7c**: GWAS overlap: bifurcation TFs and target genes × IPF GWAS hits
- **7d**: LINCS L1000 / CMap: top compounds reversing the IPF-fate signature (bar chart)
- **7e**: Proposed therapeutic model: intervention points at the bifurcation
  - Block IPF-fate TFs → prevent basaloid commitment
  - Promote repair signals → push toward AT1 maturation
  - Antagonize SPP1+ macrophage signals → remove pro-fibrotic cues

### Extended Data (≥10 supplementary figures)
- ED1: Full QC metrics, batch correction validation (kBET, LISI, silhouette)
- ED2: Complete cNMF program catalog
- ED3: Full Decoupler TF activity results (all cell types)
- ED4: SCENIC regulon analysis (full results)
- ED5: Complete LIANA+ output (all cell type pairs)
- ED6: Additional spatial validations (Mayr Xenium, Delorey GeoMx)
- ED7: Sensitivity analyses (leave-one-dataset-out)
- ED8: Donor demographics table
- ED9: Gene set enrichment details
- ED10: Code availability and reproducibility documentation

---

## 5. MANUSCRIPT NARRATIVE ARC

### Title Options (ranked)
1. **"A transcription factor switch at the KRT8+ epithelial bifurcation determines divergent repair outcomes in acute and chronic human lung injury"**
2. "Convergent injury, divergent fate: single-cell atlas analysis reveals the macrophage-epithelial decision point in COVID-ARDS versus IPF"
3. "Cross-disease atlas analysis identifies shared aberrant repair modules and a druggable fate-determining switch in human lung injury"

### Abstract Structure (~250 words)
```
BACKGROUND:  Acute (COVID-ARDS) and chronic (IPF) lung injury both destroy AT1 
             cells and activate AT2 progenitors, yet clinical outcomes diverge.
GAP:         Whether shared transitional cell states undergo identical or 
             divergent fate decisions — and what molecular switches determine 
             the outcome — remains unknown.
APPROACH:    We leveraged the harmonized HLCA (>N cells, 3 disease states) 
             to perform cross-disease trajectory, TF, and communication analysis.
FINDING 1:   KRT8+ transitional cells exist in both diseases but diverge: 
             toward apoptosis in ARDS, toward aberrant basaloid in IPF.
FINDING 2:   [Master TF] is activated at the bifurcation point exclusively 
             in IPF, locking cells into the fibrotic fate.
FINDING 3:   SPP1+ macrophages deliver divergent signals (TGFβ/WNT in IPF vs. 
             TNF/IFNγ in ARDS) that drive the fate switch.
FINDING 4:   Spatial transcriptomics confirms co-localization. The IPF-fate 
             signature predicts patient survival.
SIGNIFICANCE: Identifies a druggable fate-determining mechanism shared across 
              lung injury types with therapeutic and biomarker implications.
```

### Results Sections (7 sections matching 7 figures)
1. Harmonized atlas reveals shared and divergent cellular landscapes
2. KRT8+ transitional cells undergo fate bifurcation: apoptosis vs. basaloid
3. [Master TF] switch at the bifurcation locks IPF epithelial cells into fibrotic fate
4. Unbiased meta-programs confirm shared injury substrate with divergent overlays
5. SPP1+ macrophages deliver disease-specific fate-determining signals
6. Spatial transcriptomics validates the profibrotic niche architecture
7. The IPF-fate signature predicts survival and nominates drug targets

### Discussion (6 paragraphs)
1. Summary: first harmonized cross-disease bifurcation analysis at single-cell resolution
2. The "attractor model": why KRT8+ is a transient state in ARDS but a trap in IPF
3. Macrophage-mediated fate specification: SPP1+ cells as the "instructors"
4. Clinical implications: post-COVID fibrosis risk stratification; anti-fibrotic timing
5. Limitations (sample source bias, computational-only, no early-stage disease)
6. Future directions: experimental validation of TF KO, prospective clinical studies

---

## 6. EXECUTION: CLAUDE CODE SESSION PLAN

### Session 1 (Day 1–2): Data Extraction & QC
```
Prompt to Claude Code:
"Act as an expert bioinformatician. I am working in [Claude Code / Biomni Jupyter].
Write a Python script using cellxgene_census and scanpy to:
1. Query the Human Lung Cell Atlas
2. Filter to human cells where disease is 'normal', 'COVID-19', or 
   'idiopathic pulmonary fibrosis'
3. Filter to Epithelial + Myeloid + Stromal lineages
4. Download as AnnData, save to disk
5. Generate basic QC plots (violin: nGenes, nCounts, mito%)
6. Print cell counts per disease and cell type"
```

### Session 2 (Day 3–4): Embedding & UMAP
```
"Using the saved .h5ad from Session 1:
1. Check if scVI/scANVI latent space exists in obsm
2. If yes: compute UMAP from it; if no: run scVI with batch_key='dataset_id'
3. Generate publication-quality UMAP: 3-panel split by disease
4. Color by Level 2 cell type
5. Format for Nature: Arial 8pt, no spines, 300 DPI PDF, Set2 palette"
```

### Session 3 (Day 5–7): MiloR Compositional Analysis
```
"Using pertpy, run MiloR differential abundance:
1. Three contrasts: COVID vs. Healthy, IPF vs. Healthy, COVID vs. IPF
2. Design formula: ~ disease + sex 
3. Generate beeswarm plots for each contrast
4. Export results table with significant neighborhoods
5. Create summary Venn diagram of shared vs. specific expansions"
```

### Session 4 (Day 8–12): CellRank Trajectory ★
```
"Subset to epithelial cells. Using CellRank 2:
1. Set healthy AT2 cells as root
2. Define terminal states: AT1, Aberrant Basaloid, Apoptotic
3. Compute fate probabilities for all cells
4. Identify bifurcation point (cells with equal probability for multiple fates)
5. Generate: fate map UMAP, bifurcation score violin by disease,
   gene dynamics along branches for [KRT8, KRT17, CASP3, SOX9, COL1A1, ISG15]
6. Export lineage driver genes"
```

### Session 5 (Day 13–16): Decoupler TF Analysis ★
```
"Using Decoupler with DoRothEA (A-C confidence):
1. Compute TF activity scores for all epithelial cells
2. Subset to KRT8+ transitional cells
3. Compare TF activities: IPF-KRT8+ vs. ARDS-KRT8+ (differential activity)
4. Plot along pseudotime, split by CellRank fate
5. Generate heatmap: top 30 TFs × cell states
6. Identify the master switch TF(s)"
```

### Session 6 (Day 17–20): cNMF + LIANA+
```
"Part A — cNMF: Run consensus NMF on each compartment (Epi, Myeloid, Stromal).
Test K=5-30. Generate stability plots, program heatmaps, classification matrix.

Part B — LIANA+: Run on the full dataset with cell_type_level2 as groupby.
Focus on: SPP1+ Mac → KRT8+ Transitional, SPP1+ Mac → Fibroblast.
Compare signaling pathways across diseases."
```

### Session 7 (Day 21–26): Spatial Validation
```
"Using cell2location on Biomni GPU:
1. Train reference model on scRNA-seq (Level 2 annotations)
2. Map to IPF Visium (Franzén) and COVID Visium (Mothes)
3. Generate spatial abundance maps for key cell types
4. Run squidpy co-occurrence analysis
5. Map TF activity scores to spatial coordinates"
```

### Session 8 (Day 27–30): Clinical Validation
```
"Download GSE47460 (LGRC IPF cohort). Using lifelines:
1. Score patients for the IPF-fate TF signature
2. Generate Kaplan-Meier curves (high vs. low signature)
3. Run multivariate Cox regression
4. Query LINCS L1000 via gseapy for signature-reversing compounds
5. Generate all Figure 7 panels"
```

### Sessions 9–12 (Day 31–45): Figure Assembly & Writing
```
Iterate with Claude Code on:
- Publication-ready figure assembly (matplotlib/seaborn)
- Extended data figures
- Methods section writing
- Results section writing
- Introduction and Discussion drafts
```

---

## 7. RISK MITIGATION

| Risk | Probability | Impact | Mitigation |
|------|------------|--------|-----------|
| Census API down / rate-limited | Low | High | Fallback: download .h5ad directly from CellxGene UI |
| KRT8+ cells not identifiable in HLCA annotations | Medium | Critical | Check for: 'transitional AT2', 'KRT8+', or manually identify by marker expression (KRT8high, SFTPChigh, KRT17low) |
| CellRank can't find clean bifurcation | Medium | High | Alternative: Use Palantir (Setty 2019) for branching; or use diffusion components to visualize fate continuum; or Monocle3 |
| scVI latent space not in HLCA download | Medium | Medium | Run scVI de novo on Biomni GPU (~2–4 hours for 500K cells) |
| Survival analysis non-significant | Low-Med | Medium | Use continuous signature score (Cox) instead of binary split; try alternative cohorts (GSE70866); expand signature to include top cNMF programs |
| Memory overflow on Biomni | Medium | Medium | Subsample to 200K cells for heavy analyses (SCENIC, cNMF); use all cells for final statistics |
| Reviewer: "No experimental validation" | High | Medium | Pre-emptive: (1) Spatial validation, (2) Survival analysis, (3) GWAS overlap, (4) Cite Valenzi 2023 snATAC for chromatin-level support, (5) Offer in silico KO as "prediction for experimental follow-up" |

---

## 8. WHAT MAKES THIS v2 PLAN NATURE MEDICINE–TIER

1. **A single memorable concept**: "The KRT8+ bifurcation switch" — one sentence that editors remember
2. **Mechanistic depth**: Not just "what's different" but "which TF, activated by which macrophage signal, at which precise cell state, forces the fibrotic fate"
3. **The full causal chain**: Macrophage signal → TF activation → Fate commitment → Spatial niche → Patient survival
4. **Methodological rigor**: HLCA harmonization + orthogonal methods (Decoupler + SCENIC, LIANA+ consensus, MiloR + scCODA)
5. **Clinical hook**: Survival curves make it translational; drug candidates make it actionable
6. **The narrative arc from shared to divergent**: Convergent injury substrate → Divergent fate decision → Back to convergent implication (both diseases need fate-redirecting therapy)
7. **Visual drama**: CellRank bifurcation map with TF switches labeled is inherently a compelling figure — editors and reviewers will remember it

---

## 9. FIRST STEP TODAY

Open Claude Code or Biomni. Run Session 1. The single most important output:

> **Can you find KRT8+ transitional epithelial cells in both COVID-ARDS and IPF within the HLCA?**

If yes → the paper proceeds exactly as planned.
If the HLCA doesn't annotate them explicitly → you identify them by marker expression (KRT8+/SFTPC+/KRT17−/AGER−), which is a 20-minute scanpy operation.

Either way, you'll know within Day 1 whether this plan is viable.
