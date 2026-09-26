# Supplementary Figure 15b: scRNA-seq composition of matched cell classes

## Scripts

### `01_plot_scrna_cell_composition.R`

Reads the scRNA-seq composition table for the cell classes matched to the flow-cytometry analysis. It orders samples by tissue and age, draws the Supplementary Figure 15b stacked bars, and exports the table and PDF/PNG figure.

## Data

**[S15b_scRNAseq_Cell_Composition_source.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_15/b_S15b/03_PlotData/S15b_scRNAseq_Cell_Composition_source.csv)**

scRNA-seq counts, sample totals, and percentages for the cell classes matched to the flow-cytometry populations.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
