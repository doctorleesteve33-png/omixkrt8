# SESSION 4: CellRank Trajectory & Bifurcation Analysis
# Project: KRT8+ Epithelial Bifurcation Switch in Human Lung Injury
# Target: Nature Medicine — THIS IS FIGURE 2 (Central Figure)
# Prerequisite: Session 3 complete — MiloR + epithelial/myeloid re-clustering done

# ═══════════════════════════════════════════════════════════════════
# CONTEXT FROM SESSIONS 1-3
# ═══════════════════════════════════════════════════════════════════

"""
VERIFIED FINDINGS:
1. KRT8+ DATP cells present in BOTH COVID (~3.2% AT2) and IPF (~3.8% AT2)
2. Three disease-associated epithelial populations:
   - KRT8+ DATP: shared injury substrate (both diseases)
   - KRT5⁻KRT17+ Basaloid: IPF-specific, EMT/TGF-β/COL1A1
   - ECM-high Epithelial: COVID-specific, RBMS3/DOCK4
3. MiloR confirmed: IPF = epithelial remodeling; COVID = immune influx
4. Epithelial subset: 86,959 cells (59K Normal, 25K IPF, 2.4K COVID)
5. Key clusters: C10 (42% COVID, DATP candidate), C3 (80% IPF, basaloid candidate)
6. SPP1+ macrophages: Cluster 0, 79.5% IPF

CRITICAL CONSTRAINTS:
- Do NOT use scVelo/VelocityKernel — autopsy data has degraded spliced/unspliced
- Use DPT → PseudotimeKernel instead
- Only 2,386 COVID epithelial cells — do NOT subsample COVID further
- Harmony embedding key is 'X_harmony' (not X_pca_harmony)
- var_names must be set from 'feature_name' column before marker analysis

WHAT THIS SESSION MUST DELIVER (Figure 2):
- CellRank fate map showing KRT8+ DATP bifurcation
- Fate probabilities: P(Basaloid) vs P(ECM-high) vs P(AT1) per cell
- Bifurcation score along pseudotime
- Gene dynamics at the branching point
- Lineage driver genes for each terminal fate
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 0: DOWNLOAD + LOAD DATA
# ═══════════════════════════════════════════════════════════════════

"""
Download the 200K h5ad from GitHub Releases (same as Session 3).
Then re-derive the epithelial subset — this avoids file transfer issues.

REPO="https://github.com/doctorleesteve33-png/omixkrt8/releases/download/v1.0-data"
wget ${REPO}/hlca_part_aa through hlca_part_ah
cat hlca_part_a* > hlca_integrated_200k_subsample.h5ad
rm hlca_part_a*

Load with h5py (NOT sc.read_h5ad — causes OOM):
Use the same h5py-based loading script from Session 3.
After loading, immediately subset to epithelial cells to free memory.
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 1: REBUILD EPITHELIAL SUBSET
# ═══════════════════════════════════════════════════════════════════

"""
TASK: Subset to epithelial cells and re-process for CellRank.

  # Subset to epithelial compartment
  epi_types = ['Airway Epithelium', 'Proximal Epithelium', 
               'Alveolar Epithelium', 'Distal Epithelium',
               'AT2', 'AT1', 'Submucosal Gland', 
               'Proliferating Epi', 'Other Epi']
  # ADJUST based on actual cell_type_broad values in the data
  
  adata_epi = adata[adata.obs['cell_type_broad'].isin(epi_types)].copy()
  del adata; gc.collect()  # Free the full dataset immediately
  
  print(f"Epithelial cells: {adata_epi.shape}")
  print(adata_epi.obs['injury_type'].value_counts())
  print(adata_epi.obs['cell_type_broad'].value_counts())
  
  # Fix var_names
  if 'feature_name' in adata_epi.var.columns:
      adata_epi.var_names = adata_epi.var['feature_name'].values
      adata_epi.var_names_make_unique()
  
  # Verify key markers are present
  markers = ['KRT8', 'SFTPC', 'AGER', 'HOPX', 'KRT17', 'KRT5',
             'TP63', 'COL1A1', 'CLDN4', 'SOX4', 'CDKN1A', 'ISG15',
             'RBMS3', 'DOCK4', 'MMP7', 'LAMB3']
  found = [g for g in markers if g in adata_epi.var_names]
  print(f"Markers found: {len(found)}/{len(markers)}: {found}")

THEN RE-PROCESS:
  # HVGs on epithelial subset
  sc.pp.highly_variable_genes(adata_epi, n_top_genes=3000, 
      batch_key='source_batch', flavor='seurat_v3')
  
  # PCA
  sc.tl.pca(adata_epi, n_comps=50, use_highly_variable=True)
  
  # Harmony batch correction within epithelial subset
  import harmonypy
  ho = harmonypy.run_harmony(adata_epi.obsm['X_pca'], 
      adata_epi.obs, 'source_batch', max_iter_harmony=20)
  adata_epi.obsm['X_pca_harmony'] = ho.Z_corr.T
  
  # Neighbors + UMAP
  sc.pp.neighbors(adata_epi, use_rep='X_pca_harmony', n_neighbors=15)
  sc.tl.umap(adata_epi, min_dist=0.3)
  
  # Leiden clustering
  sc.tl.leiden(adata_epi, resolution=1.0, key_added='leiden_epi')
  
  gc.collect()
  print("Epithelial re-processing complete")
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 2: DIFFUSION PSEUDOTIME
# ═══════════════════════════════════════════════════════════════════

"""
TASK: Compute diffusion pseudotime with AT2 root cell.

WHY DPT (not scVelo): 
  Melms COVID autopsy data has degraded spliced/unspliced ratios.
  RNA velocity would produce artifacts. DPT uses expression similarity
  in diffusion space — robust to this issue.

APPROACH:
  # Diffusion map
  sc.tl.diffmap(adata_epi, n_comps=15)
  
  # Select root cell: AT2 cell in Normal tissue with highest SFTPC
  at2_normal_mask = (
      (adata_epi.obs['cell_type_broad'] == 'AT2') & 
      (adata_epi.obs['injury_type'] == 'Normal')
  )
  
  if at2_normal_mask.sum() > 0:
      # Get SFTPC expression for AT2 Normal cells
      sftpc_idx = list(adata_epi.var_names).index('SFTPC')
      from scipy.sparse import issparse
      if issparse(adata_epi.X):
          sftpc_expr = adata_epi[at2_normal_mask].X[:, sftpc_idx].toarray().flatten()
      else:
          sftpc_expr = adata_epi[at2_normal_mask].X[:, sftpc_idx].flatten()
      
      # Root = AT2 Normal cell with highest SFTPC (most "stem-like")
      at2_normal_indices = np.where(at2_normal_mask)[0]
      root_idx = at2_normal_indices[np.argmax(sftpc_expr)]
      adata_epi.uns['iroot'] = root_idx
      print(f"Root cell: index {root_idx}, SFTPC={sftpc_expr.max():.2f}")
  else:
      # Fallback: any AT2 cell with highest SFTPC
      print("WARNING: No Normal AT2 cells found, using any AT2")
      at2_mask = adata_epi.obs['cell_type_broad'] == 'AT2'
      # ... same logic
  
  # Compute DPT
  sc.tl.dpt(adata_epi)
  print(f"DPT range: {adata_epi.obs['dpt_pseudotime'].min():.3f} - "
        f"{adata_epi.obs['dpt_pseudotime'].max():.3f}")
  
  # Visualize DPT on UMAP
  sc.pl.umap(adata_epi, color=['dpt_pseudotime', 'cell_type_broad', 'injury_type'],
      save='_dpt_overview.png')

CHECKPOINT: 
  Verify DPT makes biological sense:
  - AT2 cells should have LOW pseudotime (near root)
  - AT1 cells should have MODERATE pseudotime (differentiated)
  - Basaloid/ECM-high should have HIGH pseudotime (terminal)
  If DPT looks wrong, try different root cell or n_comps.
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 3: CellRank PseudotimeKernel + GPCCA
# ═══════════════════════════════════════════════════════════════════

"""
TASK: Build CellRank transition matrix and identify terminal states.

  pip install cellrank

  import cellrank as cr
  
  # Build PseudotimeKernel from DPT
  pk = cr.kernels.PseudotimeKernel(adata_epi, time_key='dpt_pseudotime')
  pk.compute_transition_matrix(soft_threshold_scheme='exponential')
  print("Transition matrix computed")
  
  # GPCCA estimator for macrostate identification
  g = cr.estimators.GPCCA(pk)
  g.fit(cluster_key='leiden_epi', n_states=15)
  # Start with ~15 macrostates, then refine
  
  # Plot macrostates to identify terminal/initial
  g.plot_macrostates(which='all', save='cellrank_macrostates.png')
  
  # Predict terminal and initial states
  g.predict_terminal_states()
  g.predict_initial_states()
  
  print("Terminal states:", g.terminal_states.cat.categories.tolist())
  print("Initial states:", g.initial_states.cat.categories.tolist())

MANUAL VERIFICATION:
  The automatically predicted terminal states should include states 
  corresponding to:
  
  1. AT1 (repair endpoint): AGER+, HOPX+, PDPN+
  2. Aberrant Basaloid (IPF endpoint): KRT17+, COL1A1+, KRT5⁻, MMP7+
  3. ECM-Remodeling (COVID endpoint): RBMS3+, DOCK4+
  
  And the initial state should be:
  4. AT2 (progenitor): SFTPC+, ABCA3+
  
  If the automatic prediction misses any of these, manually set them:
  g.set_terminal_states(states=['state_name_1', 'state_name_2', ...])

  Check which macrostates correspond to which biology by examining:
  - Marker expression per macrostate
  - Disease composition per macrostate
  - Position on UMAP
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 4: FATE PROBABILITIES + BIFURCATION
# ═══════════════════════════════════════════════════════════════════

"""
TASK: Compute fate probabilities and bifurcation score.
THIS IS THE CORE RESULT OF THE PAPER.

  # Compute absorption probabilities to each terminal state
  g.compute_fate_probabilities()
  
  # Visualize fate probabilities on UMAP
  g.plot_fate_probabilities(save='cellrank_fate_probs.png')
  
  # ★★★ BIFURCATION SCORE ★★★
  # For each KRT8+ DATP cell: 
  # score = P(Basaloid) - P(ECM_Remodeling)
  # Positive = heading toward fibrosis
  # Negative = heading toward acute remodeling
  
  fate_probs = adata_epi.obsm['lineages_fwd']  # or g.fate_probabilities
  # Get column indices for each terminal state
  # Adjust column names based on actual terminal state names:
  basaloid_col = 'Aberrant_Basaloid'  # adjust to actual name
  ecm_col = 'ECM_Remodeling'  # adjust to actual name
  at1_col = 'AT1'  # adjust to actual name
  
  adata_epi.obs['bifurcation_score'] = (
      fate_probs[basaloid_col].values - fate_probs[ecm_col].values
  )
  
  # Plot bifurcation score on UMAP
  sc.pl.umap(adata_epi, color='bifurcation_score', 
      cmap='RdBu_r', vcenter=0, save='_bifurcation_score.png')
  
  # ★★★ BIFURCATION SCORE ALONG PSEUDOTIME ★★★
  # This is the key figure panel showing divergence
  import matplotlib.pyplot as plt
  
  fig, ax = plt.subplots(figsize=(10, 6))
  for injury, color in [('Normal', 'grey'), ('Chronic', 'red'), ('Acute', 'blue')]:
      mask = adata_epi.obs['injury_type'] == injury
      ax.scatter(
          adata_epi.obs.loc[mask, 'dpt_pseudotime'],
          adata_epi.obs.loc[mask, 'bifurcation_score'],
          c=color, alpha=0.1, s=1, label=injury
      )
  ax.axhline(y=0, color='black', linestyle='--', alpha=0.5)
  ax.set_xlabel('Pseudotime (DPT)')
  ax.set_ylabel('Bifurcation Score\nP(Basaloid) - P(ECM-Remodeling)')
  ax.legend()
  plt.savefig('bifurcation_score_pseudotime.png', dpi=300, bbox_inches='tight')
  plt.savefig('bifurcation_score_pseudotime.svg', bbox_inches='tight')
  print("Bifurcation score plot saved")

EXPECTED RESULT:
  - Early pseudotime (AT2): score near 0 (uncommitted)
  - Mid pseudotime (KRT8+ DATP): score starts diverging
  - Late pseudotime: IPF cells → positive (toward Basaloid)
  - Late pseudotime: COVID cells → negative (toward ECM-Remodeling)
  - Normal cells should cluster near 0 or toward AT1
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 5: LINEAGE DRIVERS + GENE DYNAMICS
# ═══════════════════════════════════════════════════════════════════

"""
TASK: Identify genes that drive each fate decision.

  # Compute lineage drivers (genes correlated with fate probabilities)
  drivers = g.compute_lineage_drivers(
      lineages=[basaloid_col, ecm_col, at1_col],
      return_drivers=True
  )
  
  # Save driver tables
  for lineage in [basaloid_col, ecm_col, at1_col]:
      top = drivers[lineage].sort_values(f'{lineage}_corr', ascending=False).head(50)
      top.to_csv(f'lineage_drivers_{lineage}.csv')
      print(f"\nTop 10 {lineage} drivers:")
      print(top.head(10))
  
  # ★★★ GENE DYNAMICS ALONG PSEUDOTIME ★★★
  # Plot key genes along pseudotime, split by fate
  
  key_genes = {
      'AT2_markers': ['SFTPC', 'ABCA3', 'SFTPB'],
      'DATP_markers': ['KRT8', 'CLDN4', 'SOX4', 'CEACAM6'],
      'Basaloid_markers': ['KRT17', 'COL1A1', 'MMP7', 'LAMB3'],
      'ECM_markers': ['RBMS3', 'DOCK4', 'LRMDA'],
      'AT1_markers': ['AGER', 'HOPX', 'PDPN'],
      'Arrest_markers': ['CDKN1A', 'ISG15', 'IL1B', 'GADD45A'],
      'TF_candidates': ['SOX9', 'TP63', 'HIF1A', 'STAT1', 'IRF1', 'NFKB1']
  }
  
  # For each gene category, plot expression along pseudotime
  for category, genes in key_genes.items():
      available = [g for g in genes if g in adata_epi.var_names]
      if available:
          sc.pl.scatter(adata_epi, x='dpt_pseudotime', y=available[0],
              color='injury_type', save=f'_dynamics_{category}.png')

  # CellRank gene trends along lineages
  model = cr.models.GAM(adata_epi)
  cr.pl.gene_trends(
      adata_epi, model=model,
      genes=['KRT8', 'KRT17', 'RBMS3', 'AGER', 'SFTPC', 'CDKN1A'],
      save='gene_trends_lineages.png'
  )
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 6: DISEASE-SPECIFIC FATE ANALYSIS
# ═══════════════════════════════════════════════════════════════════

"""
TASK: Compare fate probabilities between COVID and IPF KRT8+ DATP cells.
This directly tests the bifurcation hypothesis.

  # Identify KRT8+ DATP cells (from Session 1 definition)
  krt8_idx = list(adata_epi.var_names).index('KRT8')
  sftpc_idx = list(adata_epi.var_names).index('SFTPC')
  ager_idx = list(adata_epi.var_names).index('AGER')
  krt17_idx = list(adata_epi.var_names).index('KRT17')
  
  from scipy.sparse import issparse
  def get_expr(gene_idx):
      if issparse(adata_epi.X):
          return adata_epi.X[:, gene_idx].toarray().flatten()
      return adata_epi.X[:, gene_idx].flatten()
  
  krt8_expr = get_expr(krt8_idx)
  at2_median_krt8 = np.median(krt8_expr[adata_epi.obs['cell_type_broad'] == 'AT2'])
  
  datp_mask = (
      (krt8_expr > at2_median_krt8) &
      (get_expr(ager_idx) < 0.5) &
      (get_expr(krt17_idx) < 0.5)
  )
  
  print(f"KRT8+ DATP cells: {datp_mask.sum()}")
  print(f"  Normal: {(datp_mask & (adata_epi.obs['injury_type']=='Normal')).sum()}")
  print(f"  IPF: {(datp_mask & (adata_epi.obs['injury_type']=='Chronic')).sum()}")
  print(f"  COVID: {(datp_mask & (adata_epi.obs['injury_type']=='Acute')).sum()}")
  
  # Compare fate probabilities in DATP cells by disease
  datp_adata = adata_epi[datp_mask].copy()
  
  # Box/violin plot: P(Basaloid) by disease
  import matplotlib.pyplot as plt
  fig, axes = plt.subplots(1, 3, figsize=(15, 5))
  
  for i, fate in enumerate([basaloid_col, ecm_col, at1_col]):
      for injury, color in [('Normal', 'grey'), ('Chronic', 'red'), ('Acute', 'blue')]:
          mask = datp_adata.obs['injury_type'] == injury
          if mask.sum() > 0:
              axes[i].violinplot(
                  fate_probs[fate][datp_mask][mask.values],
                  positions=[['Normal','Chronic','Acute'].index(injury)],
                  showmeans=True
              )
      axes[i].set_title(f'P({fate}) in KRT8+ DATP cells')
      axes[i].set_xticks([0,1,2])
      axes[i].set_xticklabels(['Normal', 'IPF', 'COVID'])
  
  plt.tight_layout()
  plt.savefig('datp_fate_by_disease.png', dpi=300, bbox_inches='tight')
  plt.savefig('datp_fate_by_disease.svg', bbox_inches='tight')
  
  # Statistical test
  from scipy.stats import mannwhitneyu
  for fate in [basaloid_col, ecm_col, at1_col]:
      ipf_vals = fate_probs[fate][datp_mask & (adata_epi.obs['injury_type']=='Chronic')]
      covid_vals = fate_probs[fate][datp_mask & (adata_epi.obs['injury_type']=='Acute')]
      if len(ipf_vals) > 0 and len(covid_vals) > 0:
          stat, pval = mannwhitneyu(ipf_vals, covid_vals)
          print(f"{fate}: IPF median={np.median(ipf_vals):.3f}, "
                f"COVID median={np.median(covid_vals):.3f}, p={pval:.2e}")
"""

# ═══════════════════════════════════════════════════════════════════
# STEP 7: SESSION 4 OUTPUTS
# ═══════════════════════════════════════════════════════════════════

"""
MANDATORY OUTPUTS:

FIGURES (300 dpi png + svg):
  - cellrank_macrostates.png — macrostate overview
  - cellrank_terminal_states.png — terminal state UMAP
  - cellrank_fate_probs.png — fate probability per terminal state
  - bifurcation_score_umap.png — bifurcation score on UMAP
  - bifurcation_score_pseudotime.png — ★ KEY FIGURE: score vs pseudotime by disease
  - datp_fate_by_disease.png — violin: fate probs in DATP cells by disease
  - gene_trends_lineages.png — gene dynamics along each lineage
  - dpt_overview.png — pseudotime on UMAP

TABLES (csv):
  - lineage_drivers_basaloid.csv — top 50 genes driving basaloid fate
  - lineage_drivers_ecm.csv — top 50 genes driving ECM-remodeling fate
  - lineage_drivers_at1.csv — top 50 genes driving AT1 repair fate
  - datp_fate_statistics.csv — fate probabilities by disease with p-values
  - terminal_state_markers.csv — markers defining each terminal state

DATA:
  - hlca_epithelial_cellrank.h5ad — full epithelial object with fate probs
    (this feeds into Session 5 Decoupler TF analysis)

REPORT:
  Markdown summary with:
  1. How many terminal states identified? Do they match biology?
  2. Bifurcation confirmed? At what pseudotime does divergence begin?
  3. Statistical comparison: fate probs in DATP by disease
  4. Top lineage drivers for each fate
  5. Any unexpected findings or flags for Session 5
"""

# ═══════════════════════════════════════════════════════════════════
# MEMORY CONSTRAINTS
# ═══════════════════════════════════════════════════════════════════

"""
MANDATORY RULES:
  - Load 200K h5ad with h5py (NOT sc.read_h5ad)
  - Immediately subset to epithelial and delete full dataset
  - Keep X sparse throughout
  - gc.collect() after every major operation
  - If CellRank causes OOM: reduce n_states in GPCCA, or
    subsample to 50K epithelial cells (but keep ALL COVID cells)
  - Save intermediate checkpoints frequently
"""

# ═══════════════════════════════════════════════════════════════════
# DOWNSTREAM CONTEXT
# ═══════════════════════════════════════════════════════════════════

"""
SESSION 5 (Next): Decoupler TF Analysis
  - Uses hlca_epithelial_cellrank.h5ad from this session
  - Identifies the MASTER TF at the bifurcation point
  - Compares TF activities: IPF-KRT8+ vs COVID-KRT8+ cells
  - Output: siRNA targets for organoid experiments

SESSION 7: LIANA+ Cell-Cell Communication
  - Uses hlca_myeloid_subset.h5ad (from Session 3)
  - Which SPP1+ macrophage ligands predict Master TF activation?
  - Output: recombinant protein cocktail for organoid experiments
"""
