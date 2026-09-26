# Figure 3 preprocessing

| Script | Processing | Output |
| --- | --- | --- |
| `01_export_cellchat_expression.py` | Aligns raw UMI counts, log-normalized atlas expression and cell metadata. Exports sparse matrices and annotation tables for R. | `02_Processed/M03_counts.mtx.gz`, `M03_log_normalized.mtx.gz`, `M03_Cell_Metadata.csv`, `M03_Genes.tsv` |
| `02_create_seurat_input.R` | Creates a Seurat RNA assay with UMI counts and log-normalized expression. | `../../00_Global_Data/Fetal_Immune_Atlas_Seurat.rds` |
| `03_prepare_nk_velocity_layers.py` | Reads NK-cell spliced and unspliced counts, selects cells present in the atlas, and adds their cell-type and sample annotations. | `00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_07/01_NK_Trajectories/Fetal_NK_Velocity.h5ad` |
| `04_prepare_nk_cellrank_input.py` | Reads NK-cell velocity estimates, the neighborhood graph and embedding coordinates. Matches cells to the atlas annotations and assigns post-conception age to the temporal `day` field. | `00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_07/01_NK_Trajectories/Fetal_NK_CellRank_Input.h5ad` |

Run scripts 1 and 2 to generate `00_Global_Data/Fetal_Immune_Atlas_Seurat.rds` before fitting the CellChat models with `../00_CellChat_Analysis/01_Code/01_run_cellchat_for_28_cell_types.R`. The Seurat object is generated locally. Network summaries and Figure 3 panels use the supplied fitted CellChat models and plotting tables.

Scripts 3 and 4 prepare the NK-cell inputs for Supplementary Figure 7. Script 3 reads `00_Global_Data/Raw_Inputs/NK/SubNK.h5ad`, which must contain spliced and unspliced count layers. Script 4 reads `00_Global_Data/Raw_Inputs/NK/Cellrank_tpNK.h5ad`, which must contain a velocity layer, a neighborhood graph and embedding coordinates. Both scripts use the cell and sample annotations in `Fetal_Immune_Atlas.h5ad`.

The Supplementary Figure 7 scripts read the resulting h5ad files to calculate temporal transitions and draw the velocity and CellRank panels.
