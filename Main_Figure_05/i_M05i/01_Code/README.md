# Figure 5i: IGH D-J usage

`01_extract_bcr_cell_metadata.py` reads cell and sample annotations from `00_Global_Data/Fetal_Immune_Atlas_BCR.h5ad` and writes `03_Source_Data/BCR_cell_metadata.csv`.

`02_calculate_igh_dj_frequencies_and_age_correlations.R` joins these cells to the IMGT junction annotations in `00_Global_Data/Shared_Inputs/Main_Figure_04_05/05_IGH_Junction_Annotations/M05i_IGH_Junction_Annotations.txt`. It selects blood cells with an assigned D gene, calculates pooled D-J pair frequencies and within-sample percentages, and tests their association with post-conception age using Spearman correlation. BH correction is applied across the tested pairs. The complete statistics, sample percentages and 27-by-6 plotting grid are saved under `03_PlotData`.

`03_plot_igh_dj_heatmap_with_bh_significance.R` draws the frequency heatmap. Symbols show the direction of correlations with BH-adjusted P <= 0.05. Figure 5i uses IMGT D/J annotations throughout this workflow.

Run scripts 01-03 in order for the calculation and figure. Script 03 can also draw the figure from the supplied plotting table. Outputs are written to `02_Figures`.
