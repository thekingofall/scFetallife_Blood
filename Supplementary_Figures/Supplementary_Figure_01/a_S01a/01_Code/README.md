# Supplementary Figure 1a: Marker expression across cell types

## Scripts

### `01_plot_cell_type_marker_expression.py`

Reads the single-cell atlas, selects the marker genes defined in the script, and draws a stacked violin plot across the annotated cell populations. It saves Supplementary Figure 1a as PDF/PNG and exports the marker order and cell census.

## Data

**[S01a_marker_order.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_01/a_S01a/03_PlotData/S01a_marker_order.csv)**

Marker genes in the order used by the stacked violin plot.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
