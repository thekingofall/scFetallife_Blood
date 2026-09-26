# Supplementary Figure 6a: CXCR5-positive and CXCR5-negative naive B cells

## Scripts

### `01_analyze_naive_b_cell_differential_expression.py`

Reads the atlas and compares naive B-cell populations with Scanpy Wilcoxon tests. Panel a shows differential expression between CXCR5-positive and CXCR5-negative naive B cells; panel b shows organ-associated expression programs in the combined naive B-cell population. It exports ranked genes, fold-change/expression matrices, and PDF/PNG heatmaps.

## Data

**[S6a_CXCR5_NaiveB_logfoldchange_matrix.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_06/a_S06a/03_PlotData/S6a_CXCR5_NaiveB_logfoldchange_matrix.csv)**

Gene log-fold changes for the CXCR5-positive and CXCR5-negative naive B-cell comparison.

**[S6a_CXCR5_NaiveB_original_filter_top100.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_06/a_S06a/03_PlotData/S6a_CXCR5_NaiveB_original_filter_top100.csv)**

Selected differential genes from the naive B-cell comparison, with Wilcoxon scores, log-fold changes, and raw/adjusted P values.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
