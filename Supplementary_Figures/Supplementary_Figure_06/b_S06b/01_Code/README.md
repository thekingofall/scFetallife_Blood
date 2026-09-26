# Supplementary Figure 6b: Organ-associated naive B-cell expression

## Scripts

### `01_analyze_naive_b_cell_differential_expression.py`

Reads the atlas and compares naive B-cell populations with Scanpy Wilcoxon tests. Panel a shows differential expression between CXCR5-positive and CXCR5-negative naive B cells; panel b shows organ-associated expression programs in the combined naive B-cell population. It exports ranked genes, fold-change/expression matrices, and PDF/PNG heatmaps.

## Data

**[S6b_All_NaiveB_mean_expression.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_06/b_S06b/03_PlotData/S6b_All_NaiveB_mean_expression.csv)**

Mean expression of the selected naive B-cell genes in liver, spleen, and PBMCs.

**[S6b_All_NaiveB_scaled_expression.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_06/b_S06b/03_PlotData/S6b_All_NaiveB_scaled_expression.csv)**

Scaled naive B-cell expression values used for the organ-comparison heatmap.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
