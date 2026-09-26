# Figure 4i: TCR clone sizes by sample

## Scripts

### `01_plot_tcr_clone_sizes_by_sample_and_cell_type.py`

Reads the TCR clonal-expansion abundance tables and draws stacked bars grouped by sample or cell type with Scirpy. It saves the Figure 4i/4j plots as PDF and PNG and exports the ordered plotting tables.

### `02_plot_tcr_clone_sizes_by_sample.R`

Reads the sample-level TCR clone-size summary and draws stacked bars for the clone-size categories in Figure 4i. It orders samples by tissue and saves the PDF/PNG figure and source table.

## Data

**[M04i_TCR_clone_size_by_sample_data.csv](../../../00_Global_Data/Shared_Inputs/Main_Figure_04/i_M04i/03_PlotData/M04i_TCR_clone_size_by_sample_data.csv)**

Counts and proportions of clone-size categories grouped by sample and tissue.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../00_Settings/00_External_Data_Dependencies.tsv).
