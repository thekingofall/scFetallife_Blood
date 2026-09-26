# Figure 1f: Cell-lineage composition by tissue

Run the numbered scripts in order.

`01_count_cells_by_tissue_and_lineage.R`

Counts cells by tissue, cell type and lineage in the atlas metadata. Input: `00_Global_Data/Shared_Inputs/Main_Figure_01/01_Cell_Annotation/M01_Cell_Annotations_and_UMAP.csv`, exported from `Fetal_Immune_Atlas.h5ad`. The calculated table is saved in `../03_PlotData`.

`02_plot_lineage_composition.R`

Plots the eight lineage proportions within PBMC, liver, thymus and spleen. The denominator is the total number of atlas cells in each tissue. It reads the table from step 1 and saves PDF and PNG files in `../02_Figures`.
