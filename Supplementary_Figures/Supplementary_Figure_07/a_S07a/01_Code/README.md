# Supplementary Figure 7a: NK-cell RNA velocity

`01_plot_nk_rna_velocity.py` reads `00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_07/01_NK_Trajectories/Fetal_NK_CellRank_Input.h5ad`. The object contains 12,575 NK cells with the fitted velocity and UMAP coordinates used in the figure. Cell identifiers, cell-type labels and ages match the deposited cell atlas.

The script draws the scVelo velocity stream plot using the saved coordinates and the three NK-cell colors. It exports PDF and PNG files to `../02_Figures` and the cell coordinates and velocities to `../03_PlotData`.
