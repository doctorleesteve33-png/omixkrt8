# FINAL CONSOLIDATED RESEARCH PLAN v4
# The KRT8+ Epithelial Bifurcation Switch in Human Lung Injury
## Cross-Disease Single-Cell Atlas Analysis: COVID-ARDS vs. IPF

**Target journal**: Nature Medicine (primary) | Cell (backup) | Science Translational Medicine (safety)

**One-sentence pitch**: "We identify the transcription factor switch at the KRT8+ epithelial bifurcation point that determines whether injured human lung resolves (ARDS) or scars permanently (IPF), validate it epigenetically and spatially, prove causality in human alveolar organoids, and demonstrate prognostic significance for patient survival."

---

# PART 1: SCIENTIFIC FRAMEWORK

## 1.1 Two-Layer Hypothesis

**Layer 1 — Convergence (shared substrate):**
In both COVID-ARDS and IPF, AT1 cells are destroyed → AT2 progenitors become trapped in a KRT8+ transitional state → SPP1+ macrophages accumulate → CTHRC1+ fibroblasts deposit collagen. These are shared "aberrant repair modules."

**Layer 2 — Divergence (the story):**
Despite the same KRT8+ starting point, cells commit to opposite fates:
- **ARDS path**: KRT8+ → Inflammatory Arrest (CDKN1A+, ISG15+, IL1B+) → potential resolution or death
- **IPF path**: KRT8+ → Aberrant Basaloid (KRT17+, TP63low, COL1A1+) → permanent scar → no resolution

**Testable question**: Which transcription factors activate/silence at the KRT8+ bifurcation? Which macrophage-derived signals force each fate? Can we reverse the fibrotic fate?

**IMPORTANT TERMINOLOGY (v3 correction):** The ARDS terminal state is called **"Inflammatory Arrest"**, NOT "Apoptotic." scRNA-seq QC filters out truly dying cells (high mito%, low UMI). What survives QC is the pre-apoptotic inflammatory arrest state characterized by:
- Cell-cycle arrest: CDKN1A (p21), CDKN2A (p16)
- Inflammatory cytokines: IL1B, CXCL8, IL6, CCL2
- IFN response: ISG15, IFITM3, MX1, OAS1
- Stress/damage: GADD45A/B, ATF3, DDIT3

## 1.2 Why This Is Nature Medicine-Tier

1. **Single memorable concept**: "The KRT8+ bifurcation switch"
2. **Mechanistic depth**: Which TF, activated by which macrophage signal, at which cell state, forces fibrotic fate
3. **Full causal chain**: Macrophage signal → TF activation → Fate commitment → Spatial niche → Patient survival → Organoid validation → Drug rescue
4. **Methodological rigor**: HLCA harmonization + orthogonal methods (Decoupler+SCENIC, LIANA+ consensus, MiloR+scCODA)
5. **Wet-lab causal proof**: Organoid bifurcation + TF knockout + drug rescue (closes the dry-lab gap)
6. **Clinical hook**: Survival curves + drug candidates + organoid drug validation

---

# PART 2: DATASETS (All Verified Open-Access)

## 2.1 Primary scRNA-seq Data Source

**CZ CELLxGENE HLCA (Human Lung Cell Atlas)**
- Collection: https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293
- Reference: Sikkema et al. 2023, Nature Medicine
- Scale: 2.4M cells, 486 individuals, harmonized with scANVI
- Access method: CellxGene Census API (streaming, memory-efficient)

**Core datasets within HLCA:**
| Dataset | GEO | Disease | Cells | Tissue |
|---------|-----|---------|-------|--------|
| Adams et al. 2020 | GSE136831 | IPF | 312K | Lung parenchyma |
| Habermann et al. 2020 | GSE135893 | IPF | 114K | Lung parenchyma |
| Melms et al. 2021 | GSE171524 | COVID-19 | 116K | Lung autopsy |
| Delorey et al. 2021 | GSE171668 | COVID-19 | 106K | Lung autopsy |
| Reyfman et al. 2019 | GSE122960 | IPF | 76K | Lung parenchyma |

## 2.2 Spatial Transcriptomics

### Visium (spot-based, 50µm — requires cell2location deconvolution)
| Dataset | Accession | Disease |
|---------|-----------|---------|
| Franzén et al. 2024 | E-MTAB-14121 | IPF |
| Mothes et al. 2023 | Zenodo 7533977 | COVID-19 |

### Xenium (subcellular, single-molecule — NO deconvolution, use label transfer)
| Dataset | Accession | Disease |
|---------|-----------|---------|
| Vannan et al. 2025 | Per publication | IPF |
| Mayr et al. 2024 | Zenodo 10015169 | IPF |

**CRITICAL (v3 correction): Visium and Xenium require DIFFERENT analytical pipelines. See Phase 7.**

## 2.3 Epigenetic Validation (snATAC-seq / MultiOme)

**PRIMARY: Kaminski Lab 2025 MultiOme**
- GEO: **GSE286182**
- Title: "Single Nuclear RNAseq and MultiOme of IPF and Control Lung Tissue"
- Contains: snRNA-seq + snATAC-seq (MultiOme) from IPF and control lungs
- Key advantage: Has EPITHELIAL chromatin accessibility data (not just fibroblasts)
- **Use this for Phase 5b chromatin validation of the Master TF switch**

**FALLBACK: Valenzi et al. 2023**
- GEO: GSE214085
- 3 IPF + 2 control lungs, snATAC-seq
- Limitation: Primarily fibroblast-focused; epithelial coverage may be limited
- Use only if GSE286182 is inaccessible

## 2.4 Bulk RNA-seq for Clinical Validation

**LGRC Cohort**
- GEO: GSE47460
- Bulk RNA-seq from IPF patients with clinical follow-up (FVC, DLCO, survival)
- Use for: Kaplan-Meier survival analysis with Bifurcation TF Signature

**Backup**: GSE70866 (Herazo-Maya et al.)

## 2.5 Non-COVID ARDS Data (Limited Availability)

No open-access non-COVID ARDS **lung tissue** scRNA-seq exists. Available resources:
- GSE151263 (Jiang 2020): PBMCs from sepsis±ARDS patients — useful for monocyte comparison only, NO epithelial cells
- Grant et al. 2021 (GSE155249): BALF from COVID patients with some non-COVID pneumonia comparison
- **Strategy**: Check CellxGene Census for non-COVID annotations in Session 1; frame limitation transparently in Discussion; use organoid experiments to prove generalizability beyond COVID

---

# PART 3: COMPUTATIONAL PIPELINE (Phases 1–7)

## Phase 1: Data Extraction & QC (Week 1–2) — Claude Code Sessions 1–2

### Session 1: Census API Data Extraction

**Prompt for Claude Code:**
```
TASK: Extract the HLCA subset for our cross-disease study.

STEP 1 — Query Census metadata first to map exact cell_type_ontology_term_id (CL IDs):
- AT2: CL:0002063
- AT1: CL:0002062
- Alveolar macrophage: CL:0000583
- Classical monocyte: CL:0000860
- Fibroblast: CL:0000057
- Note: KRT8+ transitional and aberrant basaloid may NOT have CL IDs — 
  fall back to string matching or marker-based identification

STEP 2 — List ALL unique disease labels in the Census to check for 
non-COVID acute lung injury annotations (bacterial pneumonia, sepsis, etc.)

STEP 3 — Stream cells using cellxgene_census.get_anndata() with obs_value_filter:
- Disease conditions: 'COVID-19', 'idiopathic pulmonary fibrosis', 'normal'
- Cell compartments: Epithelial, Myeloid, Stromal
- This should yield ~200K-600K cells at ~2-5 GB memory

STEP 4 — CRITICAL VERIFICATION: Can KRT8+ transitional cells be identified 
in BOTH COVID-ARDS and IPF? Check for:
- Explicit annotation as 'transitional' or 'KRT8+' in metadata
- If not annotated: identify by marker expression (KRT8+/SFTPC+/KRT17−/AGER−)

MEMORY: Use scipy.sparse.csr_matrix everywhere. del intermediates, gc.collect().
I am running this on a cloud Jupyter with strict RAM limits.
```

**Fallback if Census API is down:** Direct .h5ad download from CellxGene UI.

### Session 2: Embedding & UMAP

**Prompt for Claude Code:**
```
TASK: Generate the unified embedding and UMAP for our cross-disease atlas.

STEP 1 — Check if the HLCA .h5ad already contains scVI/scANVI latent space 
in adata.obsm. If present, use it directly.
If absent, run scVI de novo (will need GPU — run on Biomni).

STEP 2 — Generate UMAP from the latent space.

STEP 3 — Create 3-tier cell type harmonization:
- Level 0 (Compartment): Epithelial, Myeloid, Stromal
- Level 1 (Major type): AT1, AT2, KRT8_Transitional, Basaloid, Basal, Club,
  Alv_Mac, SPP1_Mac, IFN_Mac, Classical_Mono, Fibroblast, Myofibroblast, etc.
- Level 2 (Subtype): Fine-grained subtypes per compartment

STEP 4 — QC plots: UMAP colored by disease, dataset_id, cell_type_level1, 
donor_id. Batch assessment: compute kBET and iLISI scores.

STEP 5 — Save processed .h5ad with all annotations.

MEMORY: sparse matrices, del intermediates, gc.collect().
```

## Phase 2: Compositional Analysis (Week 2–3) — Claude Code Session 3

**Prompt for Claude Code:**
```
TASK: Differential abundance analysis — which cell types expand/contract 
in each disease?

TOOL: MiloR (via milopy in Python) or pertpy

THREE CONTRASTS:
1. COVID-ARDS vs. Healthy
2. IPF vs. Healthy
3. COVID-ARDS vs. IPF

DESIGN FORMULA: ~ disease + sex + age_group
(IPF patients skew older/male; COVID patients vary. 
MUST control for age and sex as covariates.)

EXPECTED RESULTS:
- Shared expansions: SPP1+ Mac, KRT8+ Transitional, CTHRC1+ Myofibroblast
- COVID-specific: IFN-stimulated Mac, Inflammatory Monocytes
- IPF-specific: KRT17+ Aberrant Basaloid, HAS1-high Fibroblasts

VALIDATION: Run scCODA (Bayesian) as independent confirmation.

OUTPUTS:
- MiloR beeswarm plots (3 contrasts)
- Composition barplots per disease
- Venn diagram of shared vs. disease-specific abundance shifts

MEMORY: sparse matrices, del intermediates, gc.collect().
```

## Phase 3: Trajectory & Bifurcation Mapping (Week 3–5) — Claude Code Session 4 ★ CENTRAL FIGURE

**Prompt for Claude Code:**
```
TASK: Map the KRT8+ bifurcation using CellRank 2.

THIS IS THE CENTRAL ANALYSIS OF THE PAPER.

STEP 1 — Subset to EPITHELIAL compartment only.

STEP 2 — Compute RNA velocity (scVelo) or use CellRank's 
PseudotimeKernel + CytoTRACEKernel as alternatives.

STEP 3 — Define terminal states:
- Root: Healthy AT2 cells
- Terminal state 1: Mature AT1 (AGER+, HOPX+) — successful repair
- Terminal state 2: Aberrant Basaloid (KRT17+, COL1A1+, TP63low) — IPF fibrotic fate
- Terminal state 3: Inflammatory Arrest (CDKN1A+, ISG15+, IL1B+) — ARDS fate
  NOTE: Do NOT call this "Apoptotic" — see terminology note.

STEP 4 — Compute absorption probabilities for each KRT8+ cell 
toward each terminal state.

STEP 5 — BIFURCATION SCORE: For each KRT8+ cell, compute
  score = P(Basaloid) - P(Inflammatory_Arrest)
  Cells where score ≈ 0 are AT the bifurcation decision point.

STEP 6 — Gene dynamics along branches:
- ARDS-fate branch genes: CDKN1A, CDKN2A, ISG15, IFITM3, IL1B, GADD45A
- IPF-fate branch genes: KRT17, COL1A1, MMP7, SOX9, TP63

STEP 7 — CellRank lineage driver genes (volcano plot):
  genes that most strongly predict IPF-fate vs. ARDS-fate

OUTPUTS FOR FIGURE 2:
- 2a: CellRank fate map (UMAP colored by terminal state probability)
- 2b: Bifurcation score distribution (COVID vs. IPF KRT8+ cells)
- 2c: Fate probability heatmap (cells ordered by pseudotime)
- 2d: Gene dynamics along divergent branches
- 2e: Lineage driver volcano plot

MEMORY: sparse matrices, del intermediates, gc.collect().
```

## Phase 4: TF Switches at the Decision Point (Week 4–6) — Claude Code Session 5 ★ MECHANISM

**Prompt for Claude Code:**
```
TASK: Identify the Master Transcription Factor switch at the 
KRT8+ bifurcation point.

PRIMARY TOOL: Decoupler (fast, uses DoRothEA/CollecTRI databases)
SECONDARY TOOL: pySCENIC (validation, run on 50K-cell subsample on Biomni GPU)

STEP 1 — Run Decoupler MLM (multivariate linear model) on epithelial cells:
  dc.run_mlm(adata_epi, net=dorothea_abc, source='source', target='target',
             weight='weight', use_raw=False)
  Use DoRothEA confidence levels A, B, C only.

STEP 2 — Compare TF activities at the KRT8+ bifurcation:
  IPF-KRT8+ cells vs. ARDS-KRT8+ cells at same transitional state
  Use pseudobulk per donor + DESeq2 with covariates (~ disease + sex + age_group)

STEP 3 — Candidate Master Switches to specifically examine:
  IPF-fate TFs: SOX9, TP63, HIF1A, YAP1/WWTR1, TWIST1
  ARDS-fate TFs: STAT1, IRF1, NFKB1/RELA, IRF7
  Identify the ONE TF with the largest differential activity.

STEP 4 — TF activity along pseudotime, split by fate branch:
  Show the exact pseudotime window where TF activities diverge.

STEP 5 — CellOracle in silico TF knockout:
  Knock out the top candidate TF → predict shift in fate probability
  (Does KO of SOX9 shift cells from Basaloid toward AT1 repair?)

STEP 6 — pySCENIC validation (run on Biomni GPU with 64GB+ RAM):
  Subsample 50K epithelial cells, run GRNboost2 + cisTarget + AUCell
  Confirm that the Decoupler-predicted master TF has a significant regulon.

CRITICAL OUTPUT: Identify THE Master TF for:
  - Organoid Experiment B (siRNA/CRISPRi target)
  - Spatial validation (distance gradient)
  - Survival signature
  - Drug repurposing

ALSO OUTPUT: Filter the top TFs for those easily targetable by 
standard siRNA (for ordering reagents immediately).

OUTPUTS FOR FIGURE 3:
- 3a: Decoupler TF activity heatmap (IPF-KRT8+ vs. ARDS-KRT8+ vs. Healthy)
- 3b: TF activity along pseudotime, split by fate branch
- 3c: CellOracle KO simulation (fate probability shift)
- 3d: Regulon network of Master TF
- 3e: pySCENIC validation (AUCell scores)

MEMORY: sparse matrices, del intermediates, gc.collect().
```

## Phase 5: cNMF Meta-Programs + Epigenetic Validation (Week 4–6) — Claude Code Sessions 5b–6

### Session 5b: Epigenetic Chromatin Validation (NEW)

**Prompt for Claude Code:**
```
TASK: Validate the Master TF switch using chromatin accessibility 
data from snATAC-seq/MultiOme.

DATASET: GSE286182 — Kaminski Lab 2025
"Single Nuclear RNAseq and MultiOme of IPF and Control Lung Tissue"
Contains snRNA-seq + snATAC-seq (10x MultiOme) from IPF and control lungs.

IF GSE286182 IS ACCESSIBLE:
  STEP 1 — Download the MultiOme data. Focus on epithelial cells, 
  especially KRT8+ transitional and Basaloid populations.
  
  STEP 2 — Using Muon or SnapATAC2:
    - Compute TF motif deviation scores (chromVAR)
    - Compare motif accessibility: IPF epithelial vs. Control epithelial
    - Specifically test: Is our Master TF's binding motif 
      differentially accessible in IPF KRT8+/Basaloid cells?
  
  STEP 3 — Show that chromatin at the Master TF binding sites is:
    - CLOSED in healthy AT2 cells
    - OPEN in IPF KRT8+ transitional / Basaloid cells
    - This upgrades the TF finding from "computational inference" 
      to "epigenetically validated"

IF GSE286182 IS NOT ACCESSIBLE (embargo/processing):
  Use Valenzi et al. 2023 (GSE214085) as fallback.
  Note: This is fibroblast-focused, so epithelial coverage is limited.
  For fibroblasts, check TWIST1 motif accessibility as proof-of-concept.

OUTPUT: Extended Data Figure — Chromatin accessibility at Master TF 
motif sites in IPF vs. Control. Reference in Figure 3 legend.

MEMORY: sparse matrices, del intermediates, gc.collect().
```

### Session 6: cNMF Meta-Programs

**Prompt for Claude Code:**
```
TASK: Discover unbiased gene expression meta-programs using consensus NMF.

Run cNMF on EACH compartment separately:
- Epithelial (K range: 5-40, 100 iterations per K)
- Myeloid (K range: 5-40, 100 iterations per K)  
- Stromal (K range: 5-40, 100 iterations per K)

STEP 1 — K selection: Use stability vs. error plots to pick optimal K.

STEP 2 — Classify each program:
- SHARED (present in both COVID and IPF)
- COVID-SPECIFIC
- IPF-SPECIFIC
- HOMEOSTATIC (present only in healthy)

STEP 3 — Extract top 50 genes per program, run GO enrichment.

STEP 4 — Cross-compartment coordination:
  At the PATIENT level, correlate program scores across compartments.
  e.g., Does "SPP1+ Mac Program" (Myeloid) correlate with 
  "Fibrotic ECM Program" (Stromal) within the same patient?

OUTPUTS FOR FIGURE 4:
- 4a: K selection plots
- 4b: Program-by-disease classification matrix
- 4c: Top gene heatmaps per program
- 4d: Cross-compartment patient-level correlation matrix
- 4e: GO enrichment for key programs

MEMORY: cNMF is parallelizable — run on Biomni if needed.
```

## Phase 6: Cell-Cell Communication (Week 5–7) — Claude Code Session 7

**Prompt for Claude Code:**
```
TASK: Identify macrophage-derived signals that drive the KRT8+ 
bifurcation toward each fate.

PRIMARY TOOL: LIANA+ (consensus CCC framework wrapping CellPhoneDB, 
CellChat, NATMI, NicheNet — reduces false positives)

STEP 1 — Run LIANA+ on the full dataset:
  Focus on these sender→receiver pairs:
  - SPP1+ Mac → KRT8+ Transitional
  - SPP1+ Mac → Aberrant Basaloid
  - SPP1+ Mac → CTHRC1+ Fibroblast
  - IFN_Mac → KRT8+ Transitional

STEP 2 — Compare signaling pathways between diseases:
  Expected:
  - IPF-dominant: TGFβ (very high), SPP1→integrin, WNT, PDGF
  - ARDS-dominant: IL1B/TNF, IFNγ, CXCL10/CXCR3

STEP 3 — NicheNet: Which macrophage ligands predict activation of 
the Master TF discovered in Phase 4?
  - In IPF: Which ligands predict SOX9/TP63 activation?
  - In ARDS: Which ligands predict STAT1/IRF1 activation?

★★★ CRITICAL OUTPUT FOR ORGANOID EXPERIMENTS ★★★
STEP 4 — Filter the top LIANA+ predictions:
  (a) ARDS-fate ligands available as recombinant human proteins
      (e.g., TNFα, IL-1β, IFNγ — these are standard reagents)
  (b) IPF-fate ligands available as recombinant human proteins
      (e.g., TGFβ1, SPP1/OPN, WNT3A/5A)
  Output a TABLE of: Ligand | Predicted Effect | Recombinant Protein Catalog#
  This table goes directly to the bench team for ordering.

OUTPUTS FOR FIGURE 5:
- 5a: LIANA+ dot plots (SPP1+ Mac → KRT8+ Epi, by disease)
- 5b: Signaling pathway information flow comparison (COVID vs. IPF)
- 5c: NicheNet ligand-TF prediction (which ligands → which TFs)
- 5d: Ligand availability table for organoid experiments

MEMORY: sparse matrices, del intermediates, gc.collect().
```

## Phase 7: Spatial Validation & Clinical Translation (Week 6–9) — Claude Code Sessions 8–9

### Session 8: Spatial Analysis (SPLIT PIPELINE — v3 correction)

**Prompt for Claude Code:**
```
TASK: Spatial validation of the bifurcation model.

★★★ CRITICAL: Visium and Xenium require DIFFERENT pipelines ★★★

=== TRACK A: VISIUM (Deconvolution-based) ===
Datasets: Franzén IPF Visium (E-MTAB-14121), Mothes COVID Visium (Zenodo 7533977)
Method: cell2location

STEP 1 — Train reference signature model from our scRNA-seq:
  cell2location.models.RegressionModel on adata_ref with 
  labels_key='cell_type_level2', batch_key='dataset_id'

STEP 2 — Apply to Visium spatial data:
  cell2location.models.Cell2location
  Estimate per-spot abundance of each cell type.
  Run on Biomni GPU (requires ~16GB VRAM).

STEP 3 — Visualize regional gradients:
  Normal parenchyma → Transitional zone → Fibrotic core
  Show SPP1+ Mac gradient, KRT8+ Transitional gradient, Basaloid gradient

=== TRACK B: XENIUM (Direct label transfer — NO deconvolution) ===
Datasets: Vannan IPF Xenium, Mayr IPF Xenium
Method: scANVI label transfer OR scanpy.tl.ingest

STEP 1 — Transfer cell type labels from scRNA-seq to Xenium cells:
  scvi.model.SCANVI.from_scvi_model() → predict(adata_xenium)
  Each individual Xenium cell gets a cell type annotation.

STEP 2 — ★ DISTANCE-DEPENDENT SIGNALING GRADIENT (v3 addition):
  Using scipy.spatial.cKDTree:
  (a) For each KRT8+ cell, compute distance (µm) to nearest SPP1+ Mac
  (b) Extract Master TF activity score for each KRT8+ cell
  (c) Scatter plot: distance vs. TF activity
  (d) Bin by distance (0-25µm, 25-50µm, 50-100µm, 100-200µm, >200µm)
  (e) Spearman correlation + p-value
  
  EXPECTED: Exponential decay — closer to macrophage = higher TF activity.
  This mathematically proves spatially-mediated signaling.

STEP 3 — Neighborhood composition analysis:
  squidpy.gr.co_occurrence() for cell type co-localization statistics

OUTPUTS FOR FIGURE 6:
- 6a-b: Visium cell2location spatial maps (IPF and COVID tissue architecture)
- 6c: Xenium single-cell resolution co-localization
- 6d: ★ Distance-decay scatter (distance to SPP1+ Mac vs. TF activity)
- 6e: Binned distance barplot with error bars
- 6f: Co-occurrence statistics

MEMORY: cell2location needs GPU. cKDTree is memory-efficient.
```

### Session 9: Clinical Validation & Drug Repurposing

**Prompt for Claude Code:**
```
TASK: Clinical translation — survival analysis and drug candidates.

=== SURVIVAL ANALYSIS ===
Dataset: LGRC cohort GSE47460 (bulk RNA-seq with clinical follow-up)

STEP 1 — Extract "Bifurcation TF Signature":
  Top 20-50 genes that distinguish IPF-fate from ARDS-fate 
  (from CellRank lineage drivers, Phase 3)

STEP 2 — Score each patient using ssGSEA on the signature.

STEP 3 — Kaplan-Meier curves:
  Split patients by median signature score (high vs. low)
  Log-rank test for significance.

STEP 4 — Multivariate Cox proportional hazards regression:
  Adjust for age, sex, FVC%, DLCO%
  Report: Hazard Ratio, 95% CI, p-value

STEP 5 — If primary analysis is non-significant:
  Try continuous Cox score, alternative cohort (GSE70866),
  or individual TF gene expression.

=== DRUG REPURPOSING ===
STEP 6 — Query LINCS L1000 via gseapy or clue.io API:
  Input: IPF-fate gene signature (upregulated genes)
  Find compounds that REVERSE this signature.
  Prioritize FDA-approved drugs.

STEP 7 — Cross-reference with DGIdb and Open Targets 
for genetic evidence supporting druggability.

★★★ CRITICAL OUTPUT FOR ORGANOID EXPERIMENT C ★★★
STEP 8 — Select top 1-3 FDA-approved drug candidates 
for organoid rescue experiment.
Output: Drug name | Mechanism | CMap score | Clinical availability

OUTPUTS FOR FIGURE 8:
- 8a: Kaplan-Meier survival curves
- 8b: Cox regression forest plot
- 8c: GWAS overlap with bifurcation genes
- 8d: Top drug candidates from LINCS (table)
- 8e: Organoid drug rescue results (from Phase 8 wet lab)
```

---

# PART 4: ORGANOID CAUSAL VALIDATION (Phase 8 — WET LAB)

**This is the game-changer that elevates the paper from "great computational study" to "undeniable translational discovery."**

## 4.1 Overview

Three organoid experiments that physically prove the three most important computational predictions:
- **Experiment A**: Macrophage ligands dictate KRT8+ fate (tests Phase 6 LIANA+ predictions)
- **Experiment B**: Master TF is the obligatory lock for fibrotic fate (tests Phase 4 Decoupler predictions)
- **Experiment C**: Computationally predicted drug reverses fibrotic fate (tests Phase 7 CMap predictions)

## 4.2 Experiment A: Recreating the Bifurcation in a Dish (Proving Sufficiency)

**Goal:** Prove that the KRT8+ fate (Inflammatory Arrest vs. Basaloid) is strictly dictated by specific macrophage-derived ligands.

**Setup:**
1. Grow 3D human Alveolospheres from primary human AT2 cells (or iPSC-derived)
2. Establish stable baseline organoids (7-14 days culture in Matrigel)
3. Induce KRT8+ transitional state using a baseline damage trigger:
   - Option A: Mild bleomycin pulse (0.5-2 µg/mL, 24h)
   - Option B: Cytomix (TNFα + IL-1β + IFNγ at low dose)
   - Option C: Withdrawal of key growth factors to trigger AT2→transitional shift

**Divergent Perturbation (Testing LIANA+ predictions from Phase 6):**
- **Condition 1 (ARDS Simulation):** Add top ARDS macrophage ligands predicted by LIANA+
  - Expected cocktail: recombinant human TNFα + IL-1β + IFNγ (concentrations from literature: 10-50 ng/mL each)
  - **Note: Exact ligands come from Phase 6 Session 7 output**
- **Condition 2 (IPF Simulation):** Add top IPF macrophage ligands predicted by LIANA+
  - Expected cocktail: recombinant human TGFβ1 + SPP1/OPN ± WNT3A (concentrations from literature: 5-20 ng/mL)
  - **Note: Exact ligands come from Phase 6 Session 7 output**
- **Condition 3 (Control):** Vehicle only

**Readout:**
- **RT-qPCR panel:** ISG15, CDKN1A, CDKN2A, IL1B (ARDS markers) vs. KRT17, COL1A1, MMP7, ACTA2 (IPF markers) vs. AGER, HOPX (AT1 markers) vs. SFTPC (AT2 marker)
- **Immunofluorescence:** SFTPC (AT2), KRT17 (basaloid), p21/CDKN1A (arrest), AGER (AT1)
- **Expected:** Condition 1 → ISG15+/CDKN2A+ (Inflammatory Arrest); Condition 2 → KRT17+/COL1A1+ basaloid cysts

## 4.3 Experiment B: Master TF Knockout (Proving Necessity)

**Goal:** Prove that the computationally predicted Master TF is the obligatory lock that forces fibrotic fate.

**Setup:**
1. Transfect alveolar organoids with:
   - siRNA targeting the Master TF (e.g., SOX9, TP63, or HIF1A — determined by Phase 4)
   - OR lentiviral CRISPRi construct targeting the Master TF promoter
   - Scrambled siRNA/non-targeting gRNA as control
2. **Note: The specific TF target comes from Phase 4 Session 5 output**

**Perturbation:**
- Treat TF-knockdown organoids with the IPF macrophage cocktail (Condition 2 from Experiment A)

**Readout:**
- **Expected:** TF-knockdown organoids FAIL to become KRT17+ basaloid cells
- Instead, they should either:
  - Stall in KRT8+ transitional state
  - OR resolve into AGER+/HOPX+ flat AT1 cells (successful repair)
- This proves the TF is a **therapeutic bottleneck**

**Controls:**
- Scrambled siRNA + IPF cocktail → should still become KRT17+ (confirms Experiment A)
- TF-knockdown + Vehicle → baseline control
- TF-knockdown + ARDS cocktail → should still show Inflammatory Arrest (TF is IPF-specific)

## 4.4 Experiment C: CMap Drug Rescue (Translational Impact)

**Goal:** Prove that a computationally predicted FDA-approved drug reverses fibrotic fate.

**Setup:**
- IPF-cocktail-treated organoids (already becoming KRT17+)
- Add the top CMap/LINCS drug candidate at pharmacologically relevant concentration
- **Note: Drug selection comes from Phase 7 Session 9 output**

**Readout:**
- Drug should reverse KRT17+ signature
- Rescue alveolar architecture (return of SFTPC+ or AGER+ cells)
- This is the "actionable clinical translation" that top journals demand

## 4.5 Organoid Experiment Timeline

| Week | Wet Lab | Dry Lab (Claude Code) |
|------|---------|----------------------|
| 1 | Thaw AT2/iPSCs, begin expanding in 3D Matrigel | Sessions 1-3 (Census extraction, QC, embedding, MiloR) |
| 2 | Organoids expanding; order recombinant proteins + siRNAs based on Phase 4/6 outputs | Sessions 4-5 (CellRank bifurcation, Decoupler TF) |
| 3 | Organoids ready; begin damage induction | Sessions 5b-6 (MultiOme epigenetic, cNMF) |
| 4 | Run Experiment A (ARDS vs. IPF cocktails) | Session 7 (LIANA+ CCC) |
| 5 | Run Experiment B (TF knockout + IPF cocktail) | Session 8 (Spatial — Visium + Xenium) |
| 6 | Run Experiment C (Drug rescue) | Session 9 (Survival analysis, CMap) |
| 7 | Collect organoid readouts (qPCR, IF, imaging) | Session 10 (Figure assembly) |
| 8 | Organoid data analysis and figure generation | Sessions 11-12 (Manuscript writing) |

---

# PART 5: FIGURE PLAN (8 Main Figures)

## Figure 1: The Cross-Disease Unified Landscape
- 1a: Study design schematic (data sources, pipeline overview)
- 1b: UMAP split by disease (Healthy, COVID, IPF)
- 1c: MiloR beeswarm plots (3 contrasts)
- 1d: Composition barplots per disease
- 1e: Venn diagram of shared vs. disease-specific abundance shifts

## Figure 2: ★ The KRT8+ Bifurcation — CENTRAL FIGURE
- 2a: CellRank fate map (UMAP colored by terminal state probability)
- 2b: Bifurcation score distribution (COVID vs. IPF KRT8+ cells)
- 2c: Fate probability heatmap (cells ordered by pseudotime)
- 2d: Gene dynamics along divergent branches
- 2e: Lineage driver volcano plot

## Figure 3: ★ The TF Switch at the Decision Point
- 3a: Decoupler TF activity heatmap (IPF-KRT8+ vs. ARDS-KRT8+)
- 3b: TF activity along pseudotime, split by fate branch
- 3c: CellOracle in silico KO simulation
- 3d: Master TF regulon network
- 3e: pySCENIC AUCell validation

## Figure 4: cNMF Meta-Programs
- 4a: K selection plots
- 4b: Program-by-disease classification matrix
- 4c: Top gene heatmaps per program
- 4d: Cross-compartment patient-level correlation
- 4e: GO enrichment for key programs

## Figure 5: Macrophage Commands & Spatial Signaling
- 5a: LIANA+ dot plots (SPP1+ Mac → KRT8+ Epi, by disease)
- 5b: Signaling pathway information flow (COVID vs. IPF)
- 5c: NicheNet ligand-TF prediction network
- 5d: Xenium distance-decay scatter (distance to SPP1+ Mac vs. TF activity)
- 5e: Binned distance barplot

## Figure 6: Spatial Tissue Architecture
- 6a-b: Visium cell2location maps (IPF fibrotic gradient, COVID tissue)
- 6c: Xenium single-cell resolution cell type map
- 6d: Spatial co-occurrence statistics
- 6e: Niche architecture comparison (IPF vs. COVID)

## Figure 7: ★★ ORGANOID CAUSAL VALIDATION — NEW WET-LAB FIGURE
- 7a: Schematic of the hAO organoid bifurcation model
- 7b: Confocal IF images: healthy organoids (SFTPC+) vs. ARDS-cocktail (p21+/ISG15+) vs. IPF-cocktail (KRT17+/COL1A1+)
- 7c: qPCR heatmap across conditions (ARDS markers vs. IPF markers vs. AT1/AT2 markers)
- 7d: Master TF siRNA rescue: IF showing return of HOPX+/AGER+ AT1 cells upon TF knockdown
- 7e: Quantification barplots (% KRT17+, % AGER+, % p21+ per condition)

## Figure 8: Clinical Translation
- 8a: Kaplan-Meier survival curves (Bifurcation TF Signature high vs. low)
- 8b: Multivariate Cox regression forest plot
- 8c: MultiOme chromatin accessibility at Master TF binding sites (GSE286182)
- 8d: LINCS/CMap drug candidates (table)
- 8e: Organoid drug rescue (IF/qPCR showing KRT17+ reversal)
- 8f: Therapeutic model (graphical abstract of intervention strategy)

## Extended Data (≥10 figures)
- ED1: QC metrics, batch correction, kBET/iLISI scores
- ED2: Complete cell type marker dotplots
- ED3: scCODA Bayesian validation of MiloR
- ED4: Complete cNMF program catalog (all compartments)
- ED5: Full Decoupler + SCENIC results
- ED6: Complete LIANA+ interaction catalog
- ED7: Additional spatial analyses (Visium + Xenium)
- ED8: MultiOme chromatin accessibility (full TF panel)
- ED9: Organoid characterization and additional conditions
- ED10: Sensitivity analyses, demographics, donor-level variability

---

# PART 6: EXECUTION STRATEGY

## 6.1 Claude Code Session Summary

| Session | Phase | Day | Task | GPU needed? |
|---------|-------|-----|------|-------------|
| 1 | Phase 1 | 1-2 | Census API extraction, verify KRT8+ cells, check non-COVID labels | No |
| 2 | Phase 1 | 3-4 | Embedding (check scVI in HLCA or run de novo), UMAP | Maybe (scVI) |
| 3 | Phase 2 | 5-7 | MiloR compositional analysis, 3 contrasts, beeswarm plots | No |
| 4 | Phase 3 | 8-12 | ★ CellRank bifurcation, fate probabilities, gene dynamics | No |
| 5 | Phase 4 | 13-16 | ★ Decoupler TF, CellOracle KO | No |
| 5b | Phase 5 | 16-17 | MultiOme GSE286182 chromatin validation | No |
| 6 | Phase 5 | 17-20 | cNMF meta-programs | Biomni (parallelization) |
| 7 | Phase 6 | 21-24 | LIANA+ CCC, NicheNet, **output organoid reagent list** | No |
| 8 | Phase 7 | 25-28 | Spatial: Visium (cell2location on GPU) + Xenium (label transfer + distance) | Biomni (c2l) |
| 9 | Phase 7 | 29-32 | Survival (GSE47460), CMap, **output drug candidates for organoids** | No |
| 10 | — | 33-36 | Figure assembly, extended data | No |
| 11-12 | — | 37-48 | Manuscript writing | No |

## 6.2 Standard Claude Code Prompt Additions

**Add to EVERY prompt for Biomni/Claude Code:**

```python
# === STANDARD HEADER FOR ALL SESSIONS ===

# Memory management
import gc
from scipy.sparse import issparse, csr_matrix

# After every major operation:
# del intermediate_variable
# gc.collect()

# Keep matrices sparse
# if not issparse(adata.X):
#     adata.X = csr_matrix(adata.X)

# For subsetting: always .copy() to prevent views holding full data
# adata_sub = adata[mask].copy()

# Covariate control in all statistical tests:
# design = "~ disease + sex + age_group"
# IPF skews older/male; COVID varies. Must control.

# CL Ontology IDs for robust Census queries:
# AT2: CL:0002063, AT1: CL:0002062
# Alveolar Mac: CL:0000583, Classical Mono: CL:0000860
# Fibroblast: CL:0000057
# KRT8+ transitional / Basaloid: likely no CL IDs — use markers

# Terminology: "Inflammatory Arrest" NOT "Apoptotic"
```

## 6.3 Key Decision Points (Go/No-Go)

| Checkpoint | When | Go condition | No-Go action |
|-----------|------|-------------|--------------|
| KRT8+ cells exist in both diseases | Session 1 | KRT8+ identifiable by annotation or markers | Identify by marker expression (20-min operation) |
| CellRank bifurcation is clear | Session 4 | Two distinct branches from KRT8+ | Try Palantir or Monocle3 as alternatives |
| Master TF is identifiable | Session 5 | ≥1 TF with strong differential activity | Broaden to top 3 candidates, test all in organoids |
| GSE286182 accessible | Session 5b | Can download MultiOme data | Fall back to Valenzi GSE214085 |
| Survival analysis significant | Session 9 | p < 0.05 in log-rank or Cox | Try continuous score, alternative cohort, individual genes |
| Organoid bifurcation works | Experiment A | Different phenotypes in ARDS vs. IPF cocktail | Adjust concentrations, timing, damage trigger |
| TF knockdown rescues fate | Experiment B | ↓KRT17+ and/or ↑AGER+ in knockdown | Test top 3 TFs; consider combinatorial knockdown |

## 6.4 Risk Mitigation

| Risk | Probability | Mitigation |
|------|------------|------------|
| Census API down | Low | Direct .h5ad download from CellxGene UI |
| KRT8+ cells not annotated | Medium | Marker-based identification: KRT8+/SFTPC+/KRT17−/AGER− |
| scVI latent space absent in HLCA | Medium | Run scVI de novo on Biomni GPU (2-4h) |
| CellRank bifurcation unclear | Low | Palantir or Monocle3 alternatives |
| No clear single Master TF | Medium | Report top 3, test all in organoids |
| GSE286182 under embargo | Medium | Valenzi GSE214085 fallback |
| Survival non-significant | Medium | Continuous Cox score, alternative cohorts |
| Memory overflow | Medium | Subsample to 200K cells for heavy analyses |
| Non-COVID ARDS data unavailable | High | Census check + Discussion framing + organoid generalizability |
| Organoid TF knockdown fails | Low-Medium | Test top 3 TFs, try CRISPRi if siRNA insufficient |

---

# PART 7: MANUSCRIPT STRUCTURE

## Title (working)
"Divergent Transcription Factor Switches at the KRT8+ Epithelial Bifurcation Determine Resolution versus Fibrosis in Human Lung Injury"

## Abstract structure (250 words)
1. Background: Acute (ARDS) and chronic (IPF) lung injury share initial alveolar damage but diverge in outcome
2. Approach: Cross-disease single-cell atlas (>Xk cells) from harmonized HLCA, spatial transcriptomics, epigenomic validation, and organoid perturbation
3. Key finding: KRT8+ transitional cells reach a bifurcation point where [Master TF], activated by disease-specific macrophage signals, locks cells into fibrotic fate
4. Validation: Spatially validated (Xenium distance gradient), epigenetically confirmed (MultiOme chromatin), causally proven (organoid knockdown), and clinically significant (survival analysis)
5. Translation: FDA-approved drug X reverses fibrotic fate in organoids

## Discussion — key paragraphs to include
1. First cross-disease computational-experimental study of epithelial fate decisions
2. KRT8+ bifurcation model: convergent substrate, divergent mechanism
3. Master TF as therapeutic bottleneck — comparison to known biology
4. Spatial proximity proof and signaling gradient
5. **Limitation: COVID as ARDS proxy** — frame honestly, note organoid evidence supports generalizability
6. **Limitation: in vitro organoid vs. in vivo** — frame as complementary to computational predictions, not definitive
7. Clinical implications: survival signature + drug repurposing pipeline
8. Future directions: in vivo validation in mouse models, clinical trial design

---

# APPENDIX: COMPLETE DATASET ACCESSIONS

| Dataset | GEO/Accession | Type | Disease | Cells/Spots | Open? |
|---------|--------------|------|---------|-------------|-------|
| HLCA (Sikkema 2023) | CellxGene collection | scRNA-seq | Multi | 2.4M | ✅ |
| Adams 2020 | GSE136831 | scRNA-seq | IPF | 312K | ✅ |
| Habermann 2020 | GSE135893 | scRNA-seq | IPF | 114K | ✅ |
| Melms 2021 | GSE171524 | snRNA-seq | COVID | 116K | ✅ |
| Delorey 2021 | GSE171668 | snRNA-seq | COVID | 106K | ✅ |
| Reyfman 2019 | GSE122960 | scRNA-seq | IPF | 76K | ✅ |
| Ren 2021 | GSE158055 | scRNA-seq | COVID | 1.46M | ✅ |
| Franzén 2024 | E-MTAB-14121 | Visium | IPF | — | ✅ |
| Mothes 2023 | Zenodo 7533977 | Visium | COVID | — | ✅ |
| Vannan 2025 | Per publication | Xenium | IPF | 1.6M cells | ✅ |
| Mayr 2024 | Zenodo 10015169 | Xenium+Visium | IPF | — | ✅ |
| Kaminski 2025 | GSE286182 | MultiOme | IPF | — | Check |
| Valenzi 2023 | GSE214085 | snATAC-seq | IPF | — | ✅ |
| LGRC | GSE47460 | Bulk RNA-seq | IPF | — | ✅ |
| Herazo-Maya | GSE70866 | Bulk RNA-seq | IPF | — | ✅ |
| Jiang 2020 | GSE151263 | scRNA-seq | Sepsis ARDS | PBMCs only | ✅ |
| Grant 2021 | GSE155249 | scRNA-seq | COVID BALF | — | ✅ |

---

*This document supersedes v1, v2, and v3. All prior versions are archived.*
*Session 1 begins with: "Verify that KRT8+ transitional cells exist in both COVID and IPF within the HLCA Census."*
*Organoid lines should be thawed and expanded starting Week 1.*
