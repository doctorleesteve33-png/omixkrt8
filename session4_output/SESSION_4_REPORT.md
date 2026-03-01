# Session 4 Report: Trajectory & Bifurcation Analysis
## KRT8+ Epithelial Bifurcation Switch in Human Lung Injury
### Method
- Palantir trajectory analysis (Setty et al., Nature Biotechnology 2019)
- 32,386 epithelial cells (18K Normal, 12K IPF, 2.4K COVID — ALL COVID kept)
- DPT root: Normal AT2 cell with highest SFTPC expression
- 3 manually specified terminal states: AT1 (AGER+), Basaloid (KRT17+), ECM-Remodeling (RBMS3+)
### Key Results
#### 1. Bifurcation is 2-way, not 3-way
- Basaloid fate dominates (mean P=0.83), ECM and AT1 fates share identical drivers
- Effective bifurcation: Basaloid (fibrotic) vs. Alveolar/Repair (non-fibrotic)
- ECM-Remodeling did not separate as independent branch (COVID cells too few)
#### 2. Disease-specific fate bias (all epithelial cells)
- IPF cells: higher P(Basaloid), lower P(ECM) and P(AT1)
- COVID cells: lower P(Basaloid), higher P(ECM) and P(AT1)
- Bifurcation score (P_basaloid - P_ecm) clearly separates diseases on UMAP
#### 3. DATP fate comparison
- Strict DATP (KRT8-high/CLDN4+/AGER-low/SFTPC-low/KRT17-low): 1,233 cells
- COVID DATP only 11 cells — insufficient for robust comparison
- Broader KRT8+ cells (12,082): significant but small effect sizes (delta ~0.009)
- LIMITATION: COVID epithelial representation too low for cell-level fate comparison
#### 4. Lineage drivers
**Basaloid fate**: MDK (0.586), AQP5 (0.579), PRSS23 (0.540), VMO1 (0.524), SERPINB3 (0.495), LY6E (0.492), KRT19 (0.489), LCN2 (0.483)
**Alveolar/Repair fate**: PEBP4 (0.550), NAPSA (0.545), SFTPC (0.543), ABCA3 (0.532), CLDN18 (0.518), LRRK2 (0.517)
### Figures
| Figure | Description |
|--------|-------------|
| fig2_bifurcation_pseudotime | KEY: Bifurcation score vs pseudotime by disease |
| fig2_umap_panels | 4-panel: bifurcation, disease, pseudotime, entropy |
| fig2_datp_fate_violin | Fate probabilities in KRT8+ DATP by disease |
| fig2_gene_dynamics | 5-panel gene dynamics along pseudotime |
| fig2_krt8_dynamics | KRT8 expression by disease along pseudotime |
### Tables
| Table | Description |
|-------|-------------|
| palantir_fate_probs.csv | Fate probabilities for all 32K cells |
| datp_fate_statistics.csv | Broad DATP fate comparison |
| datp_fate_statistics_refined.csv | Strict DATP fate comparison |
| lineage_drivers_basaloid.csv | Basaloid lineage drivers |
| lineage_drivers_ecm.csv | ECM lineage drivers |
| lineage_drivers_at1.csv | AT1 lineage drivers |
### Flags for Session 5
1. ECM-Remodeling not separable as independent branch — use 2-way model
2. COVID epithelial cells (n=2,386) insufficient for DATP-level comparison
3. Consider supplementing with COVID BALF dataset (Grant GSE155249) in revision
4. Basaloid drivers (MDK, SERPINB3, KRT19, LCN2) → feed into Decoupler TF analysis
5. AQP5 as top basaloid driver is unexpected — may reflect proximal epithelial program
