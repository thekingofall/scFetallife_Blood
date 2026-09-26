# Figure 5e: BCR gene usage in naive B cells

## Scripts

### `01_plot_receptor_gene_usage.R`

Reads gene-usage tables and the genomic ordering of TCR/BCR gene segments. It converts counts to usage percentages, arranges genes within each segment, and draws the Figure 4e/5e heatmaps, exporting the values used in the plots.

### `02_test_gene_usage_age_and_cell_type_associations.R`

Calculates TCR/BCR gene-usage percentages, gestational-age Spearman correlations, and paired cell-subset comparisons. Age arrows use raw P < 0.05; paired Wilcoxon asterisks use BH-adjusted P < 0.05. For the BCR panel, it prepares the naive B-cell display and segment-specific scales. It exports both raw usage values and display values with the statistical annotations.

## Data

**[M05e_CXCR5pos_vs_CXCR5neg_paired_wilcoxon_BH_statistics.csv](../../../00_Global_Data/Shared_Inputs/Main_Figure_05/e_M05e/03_PlotData/M05e_CXCR5pos_vs_CXCR5neg_paired_wilcoxon_BH_statistics.csv)**

Paired cell-subset comparisons of gene-usage percentages: sample-pair counts, group medians, paired differences, Wilcoxon P values, and BH-adjusted P values.

**[M05e_gene_usage_gestational_week_spearman_statistics.csv](../../../00_Global_Data/Shared_Inputs/Main_Figure_05/e_M05e/03_PlotData/M05e_gene_usage_gestational_week_spearman_statistics.csv)**

Gene-by-gene Spearman correlations between usage and gestational age, with sample counts, raw P values, BH-adjusted values, and direction markers.

**`M05e_gene_usage_plot_data.csv`**

Gene-usage percentages by sample. `Display_value` stores the values used for the heatmap colors.

**`M05e_gene_usage_plot_data_with_gestational_week.csv`**

Gene-usage percentages by sample, cell subset, gene segment, and gestational age. `Display_value` stores the values used for the heatmap colors.

**[M05e_segment_independent_scale_ranges_PBMC.csv](../../../00_Global_Data/Shared_Inputs/Main_Figure_05/e_M05e/03_PlotData/M05e_segment_independent_scale_ranges_PBMC.csv)**

Observed gene-usage ranges and the color-scale limits for each receptor gene segment.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../00_Settings/00_External_Data_Dependencies.tsv).
