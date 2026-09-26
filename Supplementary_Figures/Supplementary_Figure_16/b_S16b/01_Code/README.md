# Supplementary Figure 16b: Expression of Olink-associated genes

Run `../a_S16a/01_Code/01_test_gene_expression_sample_order_associations.R` from the Supplementary Figure 16 directory to calculate the expression and statistics.

`02_plot_gene_expression_age_scatterplots.R` reads the resulting tables in `../03_PlotData`, plots log2(TPM + 1) against post-conception age, and adds linear fits with confidence intervals. Gene colors follow the Olink assay palette. The displayed Spearman statistics use ordered sample positions, as in the manuscript annotations; the regression lines use the plotted ages. PDF and PNG files are saved in `../02_Figures`.

The primary expression input is `00_Global_Data/Fetal_Immune_Atlas_PHA_Bulk_TPM.csv` (16 stimulated PBMC samples). `00_calculate_bulk_tpm_from_counts.R` in the panel a code folder regenerates that TPM table from the supplied counts and gene annotations.
