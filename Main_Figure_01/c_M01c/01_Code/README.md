# Figure 1c: Single-cell atlas by cell type

## Scripts

### `01_define_umap_plot_functions.R`

Defines the UMAP plotting functions used by the main atlas panels. These functions draw rasterized cell points, apply the cell-type colors, and place numbered cluster labels at cluster centers. A second function highlights selected cells against the full atlas in gray.

### `02_define_umap_arrow_layer.R`

Defines `geom_markArrow`, a ggplot2 layer for drawing the small direction arrows on the UMAP plots. Source this file before calling the layer from another script.

### `03_plot_atlas_cell_types_and_tissues.R`

Reads the shared cell-annotation UMAP table and color dictionary, and draws Figures 1c and 1d with the supplied plotting functions. Figure 1c uses cell-type colors; Figure 1d uses tissue colors. It saves both plots as PDF and PNG.

## Data

**`M01c_celltype_UMAP_plot_data.csv.gz`**

Cell-level UMAP coordinates and annotations for the tissue or cell populations shown in this panel.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../00_Settings/00_External_Data_Dependencies.tsv).
