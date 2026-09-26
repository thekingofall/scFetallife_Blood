# Supplementary Figure 3i: TCR and BCR detection in the atlas

## Scripts

### `01_define_umap_plot_functions.R`

Defines the UMAP plotting functions for the supplementary tissue and lineage panels. The functions draw full-color or gray-background embeddings and place numbered cell-type labels at adjusted cluster centers.

### `02_plot_tissue_lineage_and_receptor_umaps.R`

Reads the cell-annotation UMAP table and draws the tissue, lineage, and TCR/BCR-detection panels for Supplementary Figures 2 and 3. The TCR/BCR panel colors cells by receptor-detection status. It exports PDF/PNG figures and the corresponding cell metadata.

## Data

**[S03i_source.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_03/i_S03i/S03i_source.csv)**

Cell-level UMAP coordinates and annotations for the tissue or cell populations shown in this panel.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
