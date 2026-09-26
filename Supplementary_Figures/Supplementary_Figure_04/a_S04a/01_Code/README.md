# Supplementary Figure 4a: Expression programs during erythropoiesis

## Scripts

### `01_plot_erythroid_marker_expression.py`

Reads the atlas, selects the erythroid populations and marker groups defined in the script, and plots their expression across erythropoiesis. It exports the cell-level expression table, per-population expression summaries, and PDF/PNG heatmap.

## Data

**`S4_cell_expression_source_data.csv.gz`**

Cell-level erythroid marker expression, with sample, tissue, cell-type, and age annotations.

**[S4_expression_summary.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_04/a_S04a/03_PlotData/S4_expression_summary.csv)**

Mean log-normalized expression and the fraction of expressing cells for each erythroid population and marker gene.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
