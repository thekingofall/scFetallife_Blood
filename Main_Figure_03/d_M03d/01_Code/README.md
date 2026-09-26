# Figure 3d: Ligand-receptor interactions with NK cells

`01_prepare_hla_e_receptor_interactions.R` reads the early- and late-stage CellChat communication tables from `Main_Figure_03/00_CellChat_Analysis/02_Data/03_Plot_Ready`. It selects the eleven interactions shown in Figure 3d for three NK-cell populations and 24 sender populations, retaining significant interactions from the complete fitted networks. The plotted probability follows CellChat's -1/log(probability) transformation. The resulting table is saved in `../03_PlotData`.

`02_plot_hla_e_receptor_interactions.R` draws the six-panel dot plot from that table. Rows identify ligand-receptor pairs, columns identify senders, color represents transformed communication probability and dot size represents the permutation P-value. PDF and PNG files are saved in `../02_Figures`.
