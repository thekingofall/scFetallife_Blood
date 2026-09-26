# Figure 4e: TCR gene usage in naive CD4 and CD8 T cells

## Scripts

### `01_plot_receptor_gene_usage.R`

Reads gene-usage tables and the genomic ordering of TCR/BCR gene segments. It converts counts to usage percentages, arranges genes within each segment, and draws the Figure 4e/5e heatmaps, exporting the values used in the plots.

### `02_test_gene_usage_age_and_cell_type_associations.R`

Calculates TCR/BCR gene-usage percentages and tests their association with gestational age across PBMC samples. Age arrows use raw Spearman P < 0.05. Paired Wilcoxon comparisons between naive CD4/CD8 T cells and between CXCR5-positive/negative naive B cells use BH-adjusted P < 0.05 for asterisks. It exports the statistics, segment-specific color ranges, plotting values, and annotated heatmaps.

### `03_plot_tcr_gene_usage_with_age_annotations.R`

Reads the compact TCR gene-usage, age-correlation, and color-range tables. It draws the naive CD4/CD8 T-cell heatmaps for TRAV, TRAJ, TRBV, TRBD, and TRBJ with separate segment scales. Arrows mark age correlations with raw P < 0.05. It saves PDF/PNG figures.

## Data

**[M04e_CD4_vs_CD8_paired_wilcoxon_BH_statistics.csv](../../../00_Global_Data/Shared_Inputs/Main_Figure_04/e_M04e/03_PlotData/M04e_CD4_vs_CD8_paired_wilcoxon_BH_statistics.csv)**

Paired cell-subset comparisons of gene-usage percentages: sample-pair counts, group medians, paired differences, Wilcoxon P values, and BH-adjusted P values.

**[M04e_gene_usage_gestational_week_spearman_statistics.csv](../../../00_Global_Data/Shared_Inputs/Main_Figure_04/e_M04e/03_PlotData/M04e_gene_usage_gestational_week_spearman_statistics.csv)**

Gene-by-gene Spearman correlations between usage and gestational age, with sample counts, raw P values, BH-adjusted values, and direction markers.

**`M04e_gene_usage_plot_data.csv`**

Gene-usage percentages by sample.

**`M04e_gene_usage_plot_data_with_gestational_week.csv`**

Gene-usage percentages by sample, cell subset, gene segment, and gestational age.

**[M04e_segment_independent_scale_ranges_PBMC.csv](../../../00_Global_Data/Shared_Inputs/Main_Figure_04/e_M04e/03_PlotData/M04e_segment_independent_scale_ranges_PBMC.csv)**

Observed gene-usage ranges and the color-scale limits for each receptor gene segment.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../00_Settings/00_External_Data_Dependencies.tsv).
