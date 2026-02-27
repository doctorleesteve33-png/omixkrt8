# SESSION 3: DIFFERENTIAL ABUNDANCE + TRAJECTORY PREP
# Project: KRT8+ Epithelial Bifurcation Switch in Human Lung Injury
# Target: Nature Medicine
# Prerequisite: Session 1 completed — 200K integrated subsample + KRT8+ verified

# ═══════════════════════════════════════════════════════════════════
# CONTEXT FROM SESSION 1 (Read before executing any code)
# ═══════════════════════════════════════════════════════════════════

"""
SESSION 1 RESULTS (VERIFIED):
- 1,173,162 cells integrated (HLCA Core + Adams ILD + Melms COVID)
- 200K stratified subsample with Harmony batch correction + UMAP
- Leiden clustering at res 0.3/0.5/1.0

KRT8+ TRANSITIONAL CELLS — GREEN LIGHT:
Three disease-associated epithelial populations identified:
  1. KRT8+ DATP: ~3-4% of AT2 in BOTH COVID and IPF (shared injury substrate)
  2. KRT5⁻KRT17+ Basaloid: IPF-specific (1.53% epi), EMT/TGF-β/COL1A1
  3. ECM-high Epithelial: COVID-specific (5.74% epi), RBMS3/DOCK4

UPDATED HYPOTHESIS:
  Layer 1 (Convergence): AT2 → KRT8+ DATP (shared in both diseases) ✓
  Layer 2 (Divergence):
    - IPF: KRT8+ DATP → KRT5⁻KRT17+ Basaloid (EMT, fibrotic remodeling) ✓
    - ARDS: KRT8+ DATP → ECM-high/Metabolic Stress state ✓
      NOTE: The COVID branch is NOT the ISG15+/CDKN1A+ "Inflammatory Arrest"
      originally hypothesized. It's an ECM-remodeling state. The ISG15+/CDKN1A+
      signature may appear along the TRAJECTORY rather than as a terminal state.
      CellRank analysis in Session 4 will clarify this.

TERMINOLOGY UPDATE:
  - Keep "KRT8+ DATP" for the shared transitional state
  - "Aberrant Basaloid" for the IPF terminal state (confirmed)
  - "ECM-Remodeling State" as working name for COVID terminal state
    (will refine after CellRank trajectory analysis)

DATA CHECKPOINT:
  The .h5ad file from Session 1 should be loaded from Google Drive.
  Expected file: hlca_integrated_200k_subsample.h5ad (~4 GB)
  Contains: Harmony-corrected PCA, UMAP, Leiden clusters, cell type annotations

FULL RESEARCH PLAN: https://github.com/doctorleesteve33-png/omixkrt8
  See FINAL_RESEARCH_PLAN_v4.md for complete pipeline context.
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 0: LOAD DATA + VERIFY
# ═══════════════════════════════════════════════════════════════════

"""
TASK: Load the Session 1 checkpoint and verify integrity.

# If file is local:
import scanpy as sc
adata = sc.read_h5ad("hlca_integrated_200k_subsample.h5ad")

# If downloading from Google Drive:
# pip install gdown
# gdown "https://drive.google.com/uc?id=YOUR_FILE_ID" -O hlca_integrated_200k_subsample.h5ad --fuzzy

VERIFY:
  print(adata.shape)  # Should be ~(200000, N_genes)
  print(adata.obs.columns.tolist())  # Check available metadata
  print(adata.obsm.keys())  # Should include 'X_pca_harmony', 'X_umap'
  print(adata.obs['injury_type'].value_counts())
  print(adata.obs['disease_group'].value_counts())

ALSO CHECK (carried over from Session 1 gap):
  # Ambient RNA correction status for autopsy datasets
  autopsy_cols = [c for c in adata.obs.columns if any(x in c.lower() 
      for x in ['ambient', 'soup', 'cellbender', 'preprocess', 'qc'])]
  print(f"Ambient-related columns: {autopsy_cols}")
  
  # If no ambient columns found, check the source_batch for Melms
  melms_mask = adata.obs['source_batch'].str.contains('Melms|melms|COVID', na=False)
  print(f"Melms/COVID cells: {melms_mask.sum()}")
  # Flag in report if ambient RNA status is UNKNOWN
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 1: MiloR DIFFERENTIAL ABUNDANCE ANALYSIS
# ═══════════════════════════════════════════════════════════════════

"""
TASK: Formally quantify which cell states are significantly enriched
or depleted in each disease condition. This is Figure 1c-d in our paper.

WHY MiloR (not simple bar charts):
  - Tests abundance changes at fine-grained NEIGHBORHOOD level
  - Controls for batch effects and biological covariates
  - Provides FDR-corrected statistical significance
  - The composition shifts in Session 1 were descriptive; MiloR makes them rigorous

INSTALLATION:
  pip install milopy pertpy

APPROACH — Use pertpy's Milo wrapper (cleaner API than milopy directly):

  import pertpy as pt
  import scanpy as sc
  import numpy as np

  # Work on the full 200K subsample
  # MiloR needs the neighborhood graph
  
  # Step 1a: Build KNN graph on Harmony-corrected PCA
  sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=30)
  
  # Step 1b: Initialize Milo
  milo = pt.tl.Milo()
  mdata = milo.load(adata)
  
  # Step 1c: Make neighborhoods
  milo.make_nhoods(mdata, prop=0.1)  # Sample 10% of cells as index
  # Check: mdata['milo'].obs should have ~20K neighborhoods
  
  # Step 1d: Count cells per sample per neighborhood
  milo.count_nhoods(mdata, sample_col='donor_id')
  # CRITICAL: 'donor_id' (or 'patient_id' or 'sample_id') must be the 
  # biological replicate column, NOT batch. Check which column is available:
  # print(adata.obs[['donor_id','patient_id','sample_id']].nunique())
  
  # Step 1e: Run DA testing
  # THREE CONTRASTS (run each separately):
  
  # Contrast 1: IPF vs Normal
  milo.da_nhoods(mdata, 
      design='~ injury_type + sex + age_group',  # Control covariates
      model_contrasts='injury_typeChronic-injury_typeNormal')
  
  # Contrast 2: COVID vs Normal
  milo.da_nhoods(mdata,
      design='~ injury_type + sex + age_group',
      model_contrasts='injury_typeAcute-injury_typeNormal')
  
  # Contrast 3: COVID vs IPF (direct comparison)
  milo.da_nhoods(mdata,
      design='~ injury_type + sex + age_group', 
      model_contrasts='injury_typeAcute-injury_typeChronic')

IMPORTANT COVARIATE CHECK:
  Before running, verify that covariates exist and are usable:
  - 'sex': should be binary (M/F). Check for NaN.
  - 'age_group' or 'age': bin into groups if continuous
  - If covariates have too many NaN or are confounded with disease,
    simplify to: design='~ injury_type'
  
  # Check confounding:
  pd.crosstab(adata.obs['injury_type'], adata.obs['sex'])
  pd.crosstab(adata.obs['injury_type'], adata.obs['source_batch'])

OUTPUT FOR STEP 1:
  1. Beeswarm plot: MiloR neighborhoods on UMAP, colored by logFC
     (positive = enriched in disease, negative = depleted)
  2. Volcano-style plot: logFC vs -log10(FDR) for each contrast
  3. Table: Significantly DA neighborhoods (FDR < 0.1) annotated 
     by majority cell type
  4. Summary: Which cell types are significantly enriched in IPF vs COVID?
  
  Save: milo_results.h5mu (MuData object) for downstream use
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 2: scCODA BAYESIAN VALIDATION (OPTIONAL BUT RECOMMENDED)
# ═══════════════════════════════════════════════════════════════════

"""
TASK: Validate MiloR findings with an independent method.

scCODA uses a Bayesian model to test compositional changes while
accounting for the compositional nature of cell type proportions
(if one type goes up, others must go down — a constraint MiloR ignores).

  import pertpy as pt
  
  sccoda = pt.tl.Sccoda()
  # Aggregate to sample-level composition
  sccoda_data = sccoda.load(adata, 
      type_key='cell_type_broad',
      sample_key='donor_id',
      condition_key='injury_type')
  
  # Run model
  sccoda.prepare(sccoda_data, formula='C(injury_type)', 
      reference_cell_type='automatic')
  sccoda.run_nuts(sccoda_data, num_samples=10000)
  
  # Extract credible effects
  sccoda.summary(sccoda_data)
  # Cell types with credible inclusion probability > 0.9 are significant

OUTPUT: Table comparing MiloR vs scCODA results. Concordance = strong evidence.
If they disagree on a cell type, flag it for careful interpretation.
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 3: EPITHELIAL-FOCUSED RE-CLUSTERING
# ═══════════════════════════════════════════════════════════════════

"""
TASK: Prepare the epithelial compartment for CellRank trajectory 
analysis in Session 4. This is the setup step.

WHY: The 200K subsample includes all compartments (immune, stromal, 
endothelial). CellRank trajectory analysis needs to focus on the 
epithelial lineage only, with finer-grained clustering to resolve 
the KRT8+ DATP → Basaloid/ECM-high bifurcation.

APPROACH:
  # Subset to epithelial cells
  epi_types = ['AT2', 'AT1', 'Airway Epithelium', 'Proximal Epithelium',
               'Alveolar Epithelium', 'Basal', 'Club', 'Ciliated',
               'Goblet', 'Ionocyte', 'Neuroendocrine']
  # ADJUST the list above based on actual cell_type_broad categories in the data
  
  adata_epi = adata[adata.obs['cell_type_broad'].isin(epi_types)].copy()
  print(f"Epithelial cells: {adata_epi.shape[0]}")
  print(adata_epi.obs['cell_type_broad'].value_counts())
  print(adata_epi.obs['injury_type'].value_counts())
  
  # Re-run PCA on epithelial subset (Harmony-correct again within epi)
  sc.pp.highly_variable_genes(adata_epi, n_top_genes=3000, batch_key='source_batch')
  sc.tl.pca(adata_epi, n_comps=50, use_highly_variable=True)
  
  # Re-run Harmony on epithelial subset
  import harmonypy
  ho = harmonypy.run_harmony(adata_epi.obsm['X_pca'], 
      adata_epi.obs, 'source_batch', max_iter_harmony=20)
  adata_epi.obsm['X_pca_harmony'] = ho.Z_corr.T
  
  # Neighbors + UMAP on epithelial subset
  sc.pp.neighbors(adata_epi, use_rep='X_pca_harmony', n_neighbors=15)
  sc.tl.umap(adata_epi, min_dist=0.3)
  
  # Fine-grained Leiden clustering
  sc.tl.leiden(adata_epi, resolution=0.5, key_added='leiden_epi_0.5')
  sc.tl.leiden(adata_epi, resolution=1.0, key_added='leiden_epi_1.0')
  
  # Annotate clusters with KRT8+ DATP, Basaloid, ECM-high labels
  # Use the marker definitions from Session 1:
  #   KRT8+ DATP: KRT8↑, SFTPC↓
  #   Basaloid: KRT17+, KRT5-, COL1A1+, MMP7+
  #   ECM-high: RBMS3+, DOCK4+, OxPhos-
  
  # Generate marker dotplot for epithelial subclusters
  epi_markers = ['KRT8', 'SFTPC', 'AGER', 'HOPX',   # AT2/AT1
                 'KRT17', 'KRT5', 'TP63', 'COL1A1',   # Basaloid
                 'MMP7', 'LAMB3', 'MDK',                # Fibrotic
                 'RBMS3', 'DOCK4', 'LRMDA',             # ECM-high
                 'CDKN1A', 'ISG15', 'IL1B', 'GADD45A', # Inflammatory arrest?
                 'CLDN4', 'SOX4', 'CEACAM6',            # DATP
                 'SCGB1A1', 'FOXJ1', 'MUC5B']           # Airway
  
  sc.pl.dotplot(adata_epi, epi_markers, groupby='leiden_epi_1.0',
      save='_epi_markers_leiden.png')

OUTPUT:
  1. Epithelial UMAP colored by cell type, disease, and Leiden clusters
  2. Marker dotplot across epithelial subclusters
  3. Save: adata_epi as 'hlca_epithelial_subset.h5ad'
  4. Cell count table: epithelial subtype × disease × Leiden cluster
  
  This file will be the direct input for CellRank in Session 4.
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 4: MYELOID-FOCUSED RE-CLUSTERING (PARALLEL)
# ═══════════════════════════════════════════════════════════════════

"""
TASK: Prepare the myeloid compartment for LIANA+ analysis in Session 7.
This is a parallel setup step — not blocking the critical path.

WHY: We need to resolve SPP1+ macrophages (profibrotic, IPF-associated)
vs inflammatory macrophages (COVID-associated) for the cell-cell 
communication analysis.

APPROACH:
  # Subset to myeloid/macrophage cells
  myeloid_types = ['Myeloid', 'Macrophages', 'Monocytes', 'DC', 
                   'Dendritic cells', 'Alveolar macrophage']
  # ADJUST based on actual labels
  
  adata_mye = adata[adata.obs['cell_type_broad'].isin(myeloid_types)].copy()
  
  # Same pipeline: HVG → PCA → Harmony → UMAP → Leiden
  # Focus markers:
  mye_markers = ['SPP1', 'MARCO', 'FABP4',              # Alveolar Mac
                 'FCN1', 'S100A8', 'S100A9', 'VCAN',     # Inflammatory Mono
                 'C1QA', 'C1QB', 'TREM2',                 # TREM2+ Mac
                 'CTHRC1', 'COL1A1', 'ACTA2',             # Contamination check
                 'CD163', 'MRC1', 'MSR1',                  # M2-like
                 'IL1B', 'TNF', 'CXCL10', 'IDO1',         # M1-like / IFN
                 'MMP9', 'MMP7', 'TIMP1']                  # ECM remodeling
  
  # Look specifically for SPP1+ cluster enriched in IPF
  # This is the key macrophage subtype driving KRT8+ → Basaloid transition

OUTPUT:
  1. Myeloid UMAP colored by subtype and disease
  2. Marker dotplot identifying SPP1+ macrophage cluster
  3. Save: 'hlca_myeloid_subset.h5ad'
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 5: SESSION 3 CHECKPOINT OUTPUTS
# ═══════════════════════════════════════════════════════════════════

"""
MANDATORY SESSION OUTPUTS (save all to output directory):

FIGURES (publication-quality, 300 dpi, both .png and .svg):
  - milo_beeswarm_ipf_vs_normal.png/svg
  - milo_beeswarm_covid_vs_normal.png/svg
  - milo_beeswarm_covid_vs_ipf.png/svg
  - milo_volcano_all_contrasts.png/svg
  - epi_umap_disease.png/svg
  - epi_umap_leiden.png/svg
  - epi_dotplot_markers.png/svg
  - myeloid_umap_disease.png/svg
  - myeloid_dotplot_markers.png/svg

TABLES (.csv):
  - milo_da_results_all_contrasts.csv
  - sccoda_results.csv (if run)
  - epi_cluster_composition.csv
  - myeloid_cluster_composition.csv

DATA (.h5ad):
  - hlca_epithelial_subset.h5ad (INPUT for Session 4 CellRank)
  - hlca_myeloid_subset.h5ad (INPUT for Session 7 LIANA+)

REPORT:
  Generate a summary report (markdown) with:
  1. MiloR key findings: which cell types are significantly DA?
  2. Does MiloR confirm the descriptive composition shifts from Session 1?
  3. Epithelial re-clustering: how many subclusters? Which are KRT8+ DATP?
  4. SPP1+ macrophage identification: present in IPF? Quantify.
  5. Ambient RNA status (if checked)
  6. Any flags or issues for Session 4
"""

# ═══════════════════════════════════════════════════════════════════
# MEMORY CONSTRAINTS (Same rules as Session 1)
# ═══════════════════════════════════════════════════════════════════

"""
MANDATORY RULES:
  import gc
  from scipy.sparse import issparse, csr_matrix
  
  - After subsetting: del variables you no longer need; gc.collect()
  - Keep matrices sparse throughout
  - The 200K subsample should fit in ~8-10 GB RAM
  - Epithelial subset will be much smaller (~50-80K cells)
  - If memory issues arise, run MiloR on a further subsample (100K)
  - Save intermediate results frequently
"""

# ═══════════════════════════════════════════════════════════════════
# DOWNSTREAM CONTEXT
# ═══════════════════════════════════════════════════════════════════

"""
SESSION 4 (Next): CellRank trajectory analysis on epithelial subset
  - Uses hlca_epithelial_subset.h5ad from this session
  - CytoTRACE 2 + DPT → PseudotimeKernel (NOT scVelo)
  - Map the KRT8+ DATP → Basaloid vs ECM-high bifurcation
  - Identify lineage drivers and bifurcation point

SESSION 5: Decoupler TF analysis
  - Identify the master TF switch at the bifurcation point
  - Output: siRNA targets for organoid experiments

PARALLEL WET LAB STATUS:
  - AT2 organoids should be expanding by now
  - Waiting for Session 5 output (Master TF) to order siRNAs
  - Waiting for Session 7 output (ligands) to order recombinant proteins
"""

# ═══════════════════════════════════════════════════════════════════
# BEGIN EXECUTION
# ═══════════════════════════════════════════════════════════════════

"""
START HERE. Execute Steps 0 through 5 in order.
Report results at each checkpoint before proceeding.

PRIORITY ORDER:
  1. Load + verify data (Step 0)
  2. MiloR differential abundance (Step 1) — main deliverable
  3. Epithelial re-clustering (Step 3) — critical for Session 4
  4. scCODA validation (Step 2) — if time permits
  5. Myeloid re-clustering (Step 4) — if time permits
"""
