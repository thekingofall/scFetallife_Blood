# Figure 2a: Cell-type composition by sample

Run the numbered scripts in order.

`01_count_cells_by_sample.R`

Counts the 39 annotated cell types in each atlas sample. Input: `00_Global_Data/Shared_Inputs/Main_Figure_01/01_Cell_Annotation/M01_Cell_Annotations_and_UMAP.csv`, exported from `Fetal_Immune_Atlas.h5ad`. The calculated table is saved in `../03_PlotData`.

`02_plot_sample_cell_composition.R`

Plots one stacked bar per sample, separated by tissue, with the cell-type order and color palette used in the figure. It reads the table from step 1 and saves PDF and PNG files in `../02_Figures`.
