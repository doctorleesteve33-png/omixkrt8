# SESSION 1: DATA EXTRACTION & QC — CLAUDE CODE EXECUTION FILE
# Project: KRT8+ Epithelial Bifurcation Switch in Human Lung Injury
# Status: v4 Master Plan + v5 Amendments — LOCKED. Transitioning to EXECUTION.

# ═══════════════════════════════════════════════════════════════
# BACKGROUND CONTEXT (Read before executing any code)
# ═══════════════════════════════════════════════════════════════

"""
PROJECT OVERVIEW:
We are building a cross-disease single-cell atlas comparing COVID-ARDS and IPF
to identify the transcription factor switch at the KRT8+ epithelial bifurcation
point that determines whether injured human lung resolves or scars permanently.

TARGET JOURNAL: Nature Medicine

CORE HYPOTHESIS (Two Layers):
Layer 1 — Convergence: In both diseases, AT2 cells are destroyed → AT2 progenitors
  become trapped in KRT8+ transitional state → SPP1+ macrophages accumulate →
  CTHRC1+ fibroblasts deposit collagen (shared aberrant repair modules).
Layer 2 — Divergence: Despite same KRT8+ starting point, opposite fates:
  - ARDS path: KRT8+ → Inflammatory Arrest (CDKN1A+, ISG15+, IL1B+)
  - IPF path: KRT8+ → Aberrant Basaloid (KRT17+, COL1A1+, TP63low)

PRIMARY DATA SOURCE:
CZ CELLxGENE Human Lung Cell Atlas (HLCA)
  - Sikkema et al. 2023, Nature Medicine
  - 2.4M cells, 486 individuals, harmonized with scANVI
  - Collection: https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293

KEY DATASETS WITHIN HLCA:
  - Adams et al. 2020 (GSE136831) — IPF, 312K cells
  - Habermann et al. 2020 (GSE135893) — IPF, 114K cells
  - Melms et al. 2021 (GSE171524) — COVID autopsy, 116K cells [AUTOPSY — flag for ambient RNA]
  - Delorey et al. 2021 (GSE171668) — COVID autopsy, 106K cells [AUTOPSY — flag for ambient RNA]
  - Reyfman et al. 2019 (GSE122960) — IPF, 76K cells

TERMINOLOGY RULE:
  The ARDS terminal state is called "Inflammatory Arrest", NOT "Apoptotic."
  scRNA-seq QC filters truly apoptotic cells. What survives QC is the 
  pre-apoptotic inflammatory arrest state (CDKN1A+, ISG15+, IL1B+, GADD45A+).

FULL COMPUTATIONAL PIPELINE (Phases 1-7, this is Session 1 of ~12):
  Phase 1: Data Extraction & QC (Sessions 1-2) ← WE ARE HERE
  Phase 2: Compositional Analysis — MiloR (Session 3)
  Phase 3: Trajectory & Bifurcation — CellRank with PseudotimeKernel (Session 4)
    ★ CRITICAL: Do NOT use scVelo/VelocityKernel — autopsy data has degraded 
      spliced/unspliced ratios. Use CytoTRACE 2 + DPT → PseudotimeKernel.
  Phase 4: TF Switches — Decoupler (primary) + CellOracle + pySCENIC (Session 5)
  Phase 5: Epigenetic validation — MultiOme GSE286182 (Session 5b) + cNMF (Session 6)
  Phase 6: Cell-Cell Communication — LIANA+ (Session 7)
  Phase 7: Spatial validation (Session 8) + Clinical survival (Session 9)
"""

# ═══════════════════════════════════════════════════════════════
# STEP 0: ENVIRONMENT SETUP
# ═══════════════════════════════════════════════════════════════

"""
TASK: Set up a mutually compatible Python environment for the full pipeline.

Generate and install from a unified requirements file to avoid version conflicts.
The scanpy/scvi-tools/cellrank/anndata ecosystem is notorious for numba/jax 
conflicts if packages are installed sequentially across sessions.

Key packages needed across ALL sessions:
  - scanpy >= 1.10
  - anndata >= 0.10
  - scvi-tools >= 1.1 (for scVI, scANVI, DestVI)
  - cellrank >= 2.0
  - decoupler >= 1.6
  - pertpy >= 0.7
  - liana >= 1.0
  - cellxgene-census (pin to specific version)
  - squidpy >= 1.3
  - gseapy >= 1.1
  - muon (for MultiOme in Session 5b)
  - matplotlib, seaborn, plotly (visualization)

IMPORTANT: 
  - Do NOT install packages one-by-one across sessions
  - Resolve all dependency conflicts NOW
  - Test imports after installation
"""

# ═══════════════════════════════════════════════════════════════
# STEP 1: QUERY CENSUS METADATA (NO DATA DOWNLOAD YET)
# ═══════════════════════════════════════════════════════════════

"""
TASK: Query the CellxGene Census metadata to understand what's available
BEFORE triggering any large data downloads.

CRITICAL SAFETY RULE: Do NOT call get_anndata() until we confirm:
  (a) The metadata shapes and expected cell counts
  (b) The memory footprint is manageable (target: 2-5 GB, max: 10 GB)
  (c) We have identified the correct filter strings

Use census_version pinning for reproducibility:
  census = cellxgene_census.open_soma(census_version="2023-12-15")
  (or use the latest available LTS version — check what's available first)

SUBSTEP 1a — List ALL unique disease labels in the HLCA collection:
  Query the obs metadata column 'disease' to get all unique values.
  We need to find: 'COVID-19', 'idiopathic pulmonary fibrosis', 'normal'
  ALSO CHECK FOR: 'pneumonia', 'bacterial pneumonia', 'sepsis', 
  'acute respiratory distress syndrome', 'influenza' — any non-COVID 
  acute lung injury labels.

SUBSTEP 1b — Map cell types using CL Ontology IDs:
  Primary target CL IDs:
    AT2: CL:0002063
    AT1: CL:0002062
    Alveolar macrophage: CL:0000583
    Classical monocyte: CL:0000860
    Fibroblast: CL:0000057
  
  Also query all unique 'cell_type' string values to check for:
    - 'transitional' or 'KRT8+' (may not exist as annotation)
    - 'aberrant basaloid' or 'basaloid'
    - 'SPP1+ macrophage' or similar
  
  NOTE: KRT8+ transitional and aberrant basaloid are disease-associated 
  states that likely do NOT have standard CL ontology IDs.

SUBSTEP 1c — Count cells per disease × compartment:
  Output a pandas DataFrame with:
    Rows: disease conditions (normal, COVID-19, IPF)
    Columns: compartments (Epithelial, Myeloid, Stromal)
    Values: cell counts
  
  This tells us the memory footprint before we download.

SUBSTEP 1d — Check dataset-level metadata:
  List unique 'dataset_id' values and map them to publications.
  Specifically flag which datasets are AUTOPSY-derived (Melms, Delorey).
"""

# ═══════════════════════════════════════════════════════════════
# STEP 2: CHECK FOR AMBIENT RNA CORRECTION ON AUTOPSY DATASETS
# ═══════════════════════════════════════════════════════════════

"""
TASK: Verify whether the HLCA applied ambient RNA correction (SoupX or 
CellBender) to the autopsy datasets before integration.

WHY THIS MATTERS: Post-mortem tissue from Melms and Delorey datasets 
suffers from cell lysis. Ambient mRNA floods the suspension, causing 
healthy epithelial cells to falsely appear to express macrophage markers 
(IL1B, TNF, etc.). This would directly confound our LIANA+ cell-cell 
communication analysis in Phase 6.

HOW TO CHECK:
  - Inspect HLCA metadata columns for fields like 'preprocessing_method',
    'ambient_correction', 'cellbender_applied', or similar
  - Check the Sikkema et al. 2023 paper's methods section (the HLCA 
    harmonization pipeline likely documents this)
  - If the HLCA metadata doesn't include this info, search for it in the 
    individual dataset papers (Melms Nature 2021, Delorey Nature 2021)

DECISION:
  - If ambient correction was applied → proceed normally
  - If NOT applied or unclear → flag this as a mandatory preprocessing 
    step for Session 2. We will need to apply CellBender or SoupX to the 
    autopsy datasets before any downstream analysis.
"""

# ═══════════════════════════════════════════════════════════════
# STEP 3: VERIFY KRT8+ TRANSITIONAL CELLS IN BOTH DISEASES
# ═══════════════════════════════════════════════════════════════

"""
THIS IS THE CRITICAL GO/NO-GO CHECKPOINT FOR THE ENTIRE PROJECT.

TASK: Prove that KRT8+ transitional epithelial cells can be identified 
in BOTH COVID-ARDS and IPF within the HLCA.

APPROACH A — Check annotations:
  Search the 'cell_type' metadata for any of:
    'transitional', 'KRT8+', 'intermediate', 'damage-associated transient'
  If found in both COVID and IPF → GREEN LIGHT, proceed.

APPROACH B — If not explicitly annotated (likely scenario):
  We need to identify KRT8+ transitional cells by marker expression.
  This requires downloading a SMALL subset of epithelial cells first.
  
  Download strategy (memory-safe):
    - Filter: disease IN ('COVID-19', 'idiopathic pulmonary fibrosis', 'normal')
    - Filter: cell_type compartment = Epithelial (use appropriate filter)
    - This should be ~50K-150K cells, manageable in memory
    
  Marker-based identification of KRT8+ transitional cells:
    KRT8 HIGH (> median of AT2 cells)
    SFTPC LOW-to-MEDIUM (still has some AT2 identity)
    KRT17 NEGATIVE or LOW (not yet basaloid)
    AGER NEGATIVE (not yet AT1)
    CLDN4 HIGH (transitional marker)
    
  For each disease, count cells matching this profile.
  
  EXPECTED:
    - Normal/Healthy: Very few or none (these are injury-induced)
    - COVID-ARDS: Present (moderate numbers)
    - IPF: Present (moderate to high numbers)
    
  If KRT8+ transitional cells are found in BOTH diseases → GREEN LIGHT
  If found in only one disease → REASSESS STRATEGY (major problem)
  If found in neither → CHECK MARKERS, try alternative definitions

OUTPUT FOR THIS STEP:
  - A table: disease × KRT8+ cell count
  - A UMAP or scatter plot showing KRT8+ cells colored by disease
  - Marker dotplot: KRT8, SFTPC, AGER, KRT17, CLDN4 across cell types
"""

# ═══════════════════════════════════════════════════════════════
# STEP 4: PLAN THE FULL DATA EXTRACTION (for Session 2)
# ═══════════════════════════════════════════════════════════════

"""
Based on the metadata from Steps 1-3, plan the full data extraction 
for Session 2. Do NOT execute the full download in this session.

OUTPUT: A summary document containing:
  1. Confirmed disease labels and filter strings
  2. Cell counts per disease × compartment (from Step 1c)
  3. Estimated memory footprint for the full download
  4. Ambient RNA correction status (from Step 2)
  5. KRT8+ verification results (from Step 3)
  6. Recommended obs_value_filter string for get_anndata()
  7. Any issues or flags for Session 2

The full extraction (get_anndata for all 3 compartments × 3 diseases)
happens in Session 2, along with embedding verification and UMAP.
"""

# ═══════════════════════════════════════════════════════════════
# MEMORY & CODING CONSTRAINTS (Apply to ALL code in this session)
# ═══════════════════════════════════════════════════════════════

"""
MANDATORY RULES FOR ALL CODE:

1. MEMORY MANAGEMENT:
   import gc
   from scipy.sparse import issparse, csr_matrix
   
   - After every major operation: del intermediate_var; gc.collect()
   - Keep matrices sparse: if not issparse(adata.X): adata.X = csr_matrix(adata.X)
   - For subsetting: adata_sub = adata[mask].copy()  # .copy() prevents views
   
2. CENSUS API SAFETY:
   - ALWAYS pin census_version for reproducibility
   - ALWAYS query obs metadata as DataFrame FIRST to check counts
   - NEVER trigger get_anndata() without confirming memory footprint
   - Use obs_value_filter strings to minimize download size

3. ERROR HANDLING:
   - Wrap Census API calls in try/except
   - If Census is down, document the error and suggest fallback 
     (direct .h5ad download from CellxGene UI)

4. OUTPUT:
   - Save all intermediate results as .csv or .h5ad
   - Generate publication-quality plots (matplotlib, 300 dpi)
   - Print clear status messages at each substep
"""

# ═══════════════════════════════════════════════════════════════
# DOWNSTREAM CONTEXT (What happens after Session 1)
# ═══════════════════════════════════════════════════════════════

"""
SESSION 2 (Days 3-4): Full data extraction, verify scVI/scANVI embedding 
  in HLCA, generate UMAP, 3-tier cell type harmonization.

SESSION 3 (Days 5-7): MiloR compositional analysis (3 contrasts: 
  COVID vs Healthy, IPF vs Healthy, COVID vs IPF).
  Design formula: ~ disease + sex + age_group (MUST control covariates).

SESSION 4 (Days 8-12): ★ CellRank bifurcation mapping.
  USE PseudotimeKernel (CytoTRACE 2 + DPT). DO NOT use scVelo/VelocityKernel 
  on autopsy data.

SESSION 5 (Days 13-16): ★ Decoupler TF analysis + CellOracle KO simulation.
  Identify the Master TF. Output siRNA targets for organoid experiments.

SESSION 7 (Days 21-24): LIANA+ cell-cell communication.
  ★★★ Output: Table of ligands available as recombinant proteins for 
  organoid experiments. This goes directly to the bench team.

PARALLEL WET LAB (Starting Week 1):
  - Thaw healthy donor hAT2 cells (NOT IPF patient cells)
  - Begin expanding in 3D Matrigel
  - Order stiffened hydrogel reagents (Matrigel-Alginate or PEG-RGD, >5 kPa)
  - Order Transwell inserts (0.4 µm pore) and THP-1 cells
  - Transfect/transduce siRNA/CRISPRi in 2D BEFORE embedding in 3D
  - Wait for Session 5 output (Master TF) before ordering specific siRNAs
  - Wait for Session 7 output (ligands) before ordering specific recombinant proteins
"""

# ═══════════════════════════════════════════════════════════════
# BEGIN EXECUTION
# ═══════════════════════════════════════════════════════════════

"""
START HERE. Execute Steps 0 through 4 in order.
Report results at each checkpoint before proceeding.
"""
