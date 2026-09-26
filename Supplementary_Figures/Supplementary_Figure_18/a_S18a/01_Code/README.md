# Supplementary Figure 18a: MOFA factor 3 feature weights

`01_plot_mofa_factor_feature_weights.R` reads [S18a_Factor3_Top6_Weights_source.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_18/a_S18a/03_PlotData/S18a_Factor3_Top6_Weights_source.csv). For each of the 25 data views, it plots the six largest positive and six largest negative normalized feature weights for factor 3. View names and colors follow the manuscript figure. PDF and PNG files are saved in `../02_Figures`.

Run the script with Rscript. Paths are resolved from the script location; R libraries are configured through `00_Settings/00_setup_R.R`.
