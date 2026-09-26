# Figure 1e: Tissue distribution of cell types

Run the numbered scripts in order.

`01_calculate_tissue_cell_type_fractions.R`

Counts each cell type per sample and calculates its sample-level fraction. Tissue summaries use a weighted mean, with weights equal to the sample cell count multiplied by the number of samples from that tissue. Samples with no observed cells of a given type are absent from that type's grouping. Input: `00_Global_Data/Shared_Inputs/Main_Figure_01/01_Cell_Annotation/M01_Cell_Annotations_and_UMAP.csv`, exported from `Fetal_Immune_Atlas.h5ad`. The calculated table is saved in `../03_PlotData`.

`02_plot_tissue_cell_type_fractions.R`

Plots the tissue fractions as stacked bars, grouped by cell lineage, using the manuscript tissue colors. It reads the table from step 1 and saves PDF and PNG files in `../02_Figures`.
