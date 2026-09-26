# Figure 4j: TCR clone sizes by cell type

## Scripts

### `01_plot_tcr_clone_sizes_by_sample_and_cell_type.py`

Reads the TCR clonal-expansion abundance tables and draws stacked bars grouped by sample or cell type with Scirpy. It saves the Figure 4i/4j plots as PDF and PNG and exports the ordered plotting tables.

### `02_plot_tcr_clone_sizes_by_cell_type.R`

Reads the cell-type-level TCR clone-size summary and plots the distribution of clone-size categories across T-cell populations for Figure 4j. It exports the plotted table and PDF/PNG figures.

## Data

**[M04j_TCR_clone_size_by_celltype_data.csv](../../../00_Global_Data/Shared_Inputs/Main_Figure_04/j_M04j/03_PlotData/M04j_TCR_clone_size_by_celltype_data.csv)**

Counts and proportions of clone-size categories grouped by cell type and tissue.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../00_Settings/00_External_Data_Dependencies.tsv).
