# Circulating immune programs and TCR repertoire associated with survival after CMV-specific ACT in recurrent glioblastoma

Analysis code for single-cell RNA-seq, CITE-seq and TCR-seq of PBMCs collected before (Pre) and after (Post) CMV-specific adoptive T-cell therapy (ACT) in 9 patients with recurrent GBM.

## Repository layout

```
R/helpers.R                     shared paths, palettes and helper functions
scripts/01_quality_control.R    per-sample QC, doublet removal, cell cycle, merge
scripts/02_TCR_preprocessing.R  VDJ contigs -> one TCR per cell, added to Seurat
scripts/03_integration.R        SCTransform + Harmony (sequencing run), contaminant removal
scripts/04_clustering.R         clusters, cell-type annotation, subclusters, DE
scripts/05_hdWGCNA.R            co-expression network, module scores, module figures
scripts/06_TCR_repertoire.R     D50, clonal expansion, pre/post sharing, ACT-derived clones
scripts/07_correlation_analysis.R  patient-level associations with overall survival
data/metadata/                  input metadata (see data/metadata/README.md)
```

Run the scripts in order from the repository root, e.g. `Rscript scripts/01_quality_control.R`. Each script reads the `.rds` written by the previous one from `results/rds/` and writes tables to `results/tables/` and figures to `results/figures/`. Set the environment variable `GBM_PROJECT_DIR` if data and results live elsewhere.

## Input data

- Cell Ranger `multi` outputs in `data/raw/GBM*/` (gene expression, antibody capture, VDJ-T). Raw data will be available from EGA after archiving process is complete.

## Figure to script map

| Figure | Script |
|---|---|
| 2A | 04_clustering.R |
| 2B, Supp 1A, 1B | 07_correlation_analysis.R |
| 2C, Supp 1D | 05_hdWGCNA.R (scores), 07_correlation_analysis.R |
| 2D, 2E, Supp 1C, Supp 2 | 05_hdWGCNA.R |
| 3A, 3C, Supp 3A | 04_clustering.R |
| 3B | 05_hdWGCNA.R |
| 4A, Supp 3B, Supp 4 | 07_correlation_analysis.R |
| 4B | 04_clustering.R (DEGs), 05_hdWGCNA.R (Venn) |
| 4C | 04_clustering.R |
| 5A, 5C, 6A, 6C, 7A, Supp 6 | 06_TCR_repertoire.R |
| 5B, 6B, 7C, Supp 5 | 06_TCR_repertoire.R (tables), 07_correlation_analysis.R |
| 7B | 04_clustering.R |

## Key definitions

- **Clonotype:** TRB V gene + J gene + CDR3 amino-acid sequence. Cells with two TRB (or TRA) chains keep the chain with the most UMIs; ties are removed.
- **D50:** percentage of unique clonotypes that together account for 50% of TRB+ cells.
- **ACT-derived clones:** single-cell TRB clonotypes matching the TCR repertoire of the patient's ACT product.
- **Module scores:** `AddModuleScore` on genes of each hdWGCNA module with kME ≥ 0.5.
- **Statistics:** Pearson correlation with OS from diagnosis (n = 9).

## Software

R 4.x with Seurat v5, harmony, DoubletFinder, hdWGCNA, WGCNA, enrichR, circlize, VennDiagram, pheatmap, tidyverse. Exact versions are printed by `sessionInfo()` at the end of `07_correlation_analysis.R`.

## Citation

