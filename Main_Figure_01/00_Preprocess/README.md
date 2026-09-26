# Figure 1 preprocessing

Run the Python scripts in numeric order from any working directory.

| Script | Processing | Output |
| --- | --- | --- |
| `01_score_doublets_and_predict_cell_types.py` | Reads the 31 sample matrices listed in `01_Inputs/Sample_Inputs.csv`. Adds the sample prefix to cell barcodes, scores doublets with Scrublet and calls doublets at a fixed threshold of 0.8 for each sample, merges the samples, normalizes to 10,000 counts per cell, applies log1p and assigns CellTypist labels. | `02_Processed/M01_Raw_Counts.h5ad`, `M01_CellTypist.h5ad` |
| `02_filter_cells_and_cluster.py` | Adds sample metadata and TCR/BCR detection labels; removes predicted doublets and cells detected in both receptor assays. Applies the recorded RNA quality filters, selects variable genes, regresses total counts and mitochondrial percentage, scales expression, and calculates PCA, neighbors, UMAP and Leiden clusters. | `02_Processed/M01_Filtered_Clustered.h5ad` |
| `03_apply_cell_annotations.py` | Matches the retained cells to the curated cell annotations, gene list and manuscript UMAP coordinates. Writes the full-gene, log-normalized atlas with cell and sample metadata. | `../../00_Global_Data/Fetal_Immune_Atlas.h5ad` |

## Inputs

Place the Cell Ranger matrices at the relative paths in `01_Inputs/Sample_Inputs.csv`. Place the CellTypist model at `../../00_Global_Data/Raw_Inputs/CellTypist/Immune_All_Low.pkl`. TCR/BCR annotations are read from `00_Global_Data/Shared_Inputs/Main_Figure_04_05/00_Receptor_Annotations/01_Aggregated_Contigs`.

The curated cell assignments and manuscript coordinates are in `00_Global_Data/Shared_Inputs/Main_Figure_01/01_Cell_Annotation/M01_Cell_Annotations_and_UMAP.csv`. `01_Inputs/Cell_Type_Markers.json` contains the marker lists used to examine the progenitor, erythroid, myeloid, T-cell and B-cell annotations. Manual cell assignments are an input to step 3.

The quality filters are: at least 200 detected genes per cell; genes detected in at least 3 cells; fewer than 8,000 detected genes per cell; mitochondrial percentage below 10. QC metrics are calculated on the normalized, log-transformed matrix. PCA uses the ARPACK solver; the neighbor graph uses 10 neighbors and 40 PCs.

The specified Python versions are in `../../00_Settings/Preprocess_Python_Versions.tsv`.

The figure scripts read `00_Global_Data/Fetal_Immune_Atlas.h5ad`. The preprocessing scripts write intermediate files to `02_Processed` and the annotated atlas to `00_Global_Data`.
