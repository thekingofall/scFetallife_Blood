# Supplementary Figure 3g: Megakaryocyte and erythroid cells

## Scripts

### `01_define_umap_plot_functions.R`

Defines the rasterized UMAP plotting functions for tissue and lineage panels. It calculates cluster centers, separates overlapping numbered labels, and supports colored subsets over the gray atlas background.

### `02_plot_tissue_lineage_and_receptor_umaps.R`

Reads the shared UMAP coordinates, cell annotations, and palette. It draws the four tissue panels in Supplementary Figure 2, the eight lineage panels in Supplementary Figure 3a-h, and the TCR/BCR-detection panel in Supplementary Figure 3i. Each panel is saved as PDF/PNG with its selected-cell table.

## Data

**[S03g_source.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_03/g_S03g/S03g_source.csv)**

Cell-level UMAP coordinates and annotations for the tissue or cell populations shown in this panel.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
