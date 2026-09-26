# Supplementary Figure 08b: CXCR6 NK-cell ligand–receptor interactions

`01_plot_nk_ligand_receptor_heatmap.R` reads [S08b_CXCR6_NK_Interaction_Matrix.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_08/03_PlotData/S08b_CXCR6_NK_Interaction_Matrix.csv) and plots the prior interaction potential for each ligand–receptor pair with NicheNet. Receptors are columns and ligands are rows. The heatmap is exported to `../02_Figures` as PDF and PNG.

Run the script with Rscript. Input and output paths are resolved from the script location. The script uses the shared R settings in `00_Settings`.
