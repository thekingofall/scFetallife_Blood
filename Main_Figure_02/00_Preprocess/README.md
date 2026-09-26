# Figure 2 preprocessing

The cell atlas is produced by `Main_Figure_01/00_Preprocess`. The scripts below prepare inputs for the gene-module analysis and the fetal/reference comparison.

| Script | Processing | Output |
| --- | --- | --- |
| `01_prepare_hotspot_expression.py` | Selects blood cells and the recorded variable-gene set; removes mitochondrial and ribosomal genes; aligns raw UMI counts, normalized expression and stored PCA coordinates. Splits samples at 26 post-conception weeks. | `02_Processed/M02_PBMC_Early_Hotspot.h5ad`, `M02_PBMC_Late_Hotspot.h5ad` |
| `02_calculate_gene_correlations.py` | Fits the Hotspot DANB model using 25 neighbors. Selects genes with autocorrelation FDR below 0.00001, computes local correlations and identifies modules with at least 20 genes and FDR 0.01. | `early_modules.csv`, `late_modules.csv`, their local-correlation matrices under `00_Global_Data/Shared_Inputs/Main_Figure_02/c_M02c`, and `M02d_Early_Detected_Modules.csv` / `M02d_Late_Detected_Modules.csv` under `d_M02d` |
| `03_prepare_reference_metadata_h5ad.py` | Combines the curated fetal/reference metadata with the saved UMAP coordinates in `01_Inputs/Reference_Cell_Metadata_and_UMAP.csv.gz`. The cell-type palette is read from `01_Inputs/Reference_Colors.json`. The resulting h5ad has 184,222 cells and no expression columns. | `00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/01_Cell_Metadata/S05_reference_metadata_umap.h5ad` |

Step 1 requires the raw sample matrices listed in the Figure 1 input manifest. Run step 2 with `--jobs 8` to set its worker count. The downstream `d_M02d` scripts map detected modules to the manuscript module names and perform enrichment.
