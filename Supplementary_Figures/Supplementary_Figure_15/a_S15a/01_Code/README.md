# Supplementary Figure 15a: spectral flow-cytometry cell composition

`01_aggregate_and_plot_flow_cell_composition.R` reads the blood and organ CD45/lymphocyte frequency tables from `00_Global_Data/Shared_Inputs/Main_Figure_06/00_Flow_Summaries`. It combines the 17 gated populations into ten measured cell classes and calculates the remaining percentage as Ungate. The panel contains 27 samples from blood, liver, thymus and spleen.

Blood percentages use the lymphocyte parent population; organ percentages use CD45+ cells. The B17.4_P3 DNT percentage is read from the unrounded flow-cytometry workbook under `Shared_Inputs/Supplementary_Figures/Supplementary_Figure_15/a_S15a`.

The script saves the 17-population input table and the aggregated composition tables in `../03_PlotData`, then exports the stacked-bar figure to `../02_Figures` as PDF and PNG. Run with Rscript; paths are resolved from the script location.
