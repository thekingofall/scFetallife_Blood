# Figure 7b: Selected plasma proteins

`01_plot_selected_plasma_proteins.R`

Plots NPX values against age for 17 selected plasma proteins across 18 samples. Points represent individual samples and lines show linear fits with confidence bands. NPX measurements are read from `00_Global_Data/Shared_Inputs/Main_Figure_07/M07b_Selected_Plasma_Protein_NPX.csv`. Spearman r, P and BH-adjusted q values are read from `M07a_Plasma_Protein_Age_Statistics.tsv` in the same directory, sharing the 76-protein statistical results with Figure 7a.

PDF and PNG files are saved in `../02_Figures`. Plotting tables are exported to `00_Global_Data/Shared_Inputs/Main_Figure_07/b_M07b/03_PlotData`.
