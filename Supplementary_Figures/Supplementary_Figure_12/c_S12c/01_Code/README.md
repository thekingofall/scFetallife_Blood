# Supplementary Figure 12c: IGL junction lengths

`01_calculate_junction_length_age_associations.R` reads the IMGT-derived sample-by-length frequency tables in `00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_12/01_Junction_Length`. It extracts post-conception age from the sample ID and applies the Spearman test to the frequency of each nucleotide length in each junction component. Distribution data and correlation statistics are saved in `../03_PlotData`.

`02_plot_junction_length_distributions_and_correlations.R` draws each junction component in one row. The left plot shows sample-level frequencies and smoothed curves colored by age. The right plot shows Spearman r with the figure's P-value categories and palette. PDF and PNG files are saved in `../02_Figures`.
