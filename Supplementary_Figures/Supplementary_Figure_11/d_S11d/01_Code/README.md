# Supplementary Figure 11d: Selected TCR clones by sample

Run the network script for panel b first. `01_count_cells_in_selected_tcr_clones.py` reads the resulting clone assignments and the TCR atlas cell types, then counts cells in the 15 clonotypes selected for the manuscript. The selected clone list is stored under `00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_11/d_S11d`. It writes the sample-by-cell-type counts to `../03_PlotData`.

`02_plot_top_tcr_clones_by_sample.R` draws the clone composition bars in seven sample facets, with gray strips, black borders and cell-type colors. The panel contains 73 cells. PDF and PNG files are saved in `../02_Figures`.
