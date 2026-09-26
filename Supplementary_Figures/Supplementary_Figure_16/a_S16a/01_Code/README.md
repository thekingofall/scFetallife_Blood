# Supplementary Figure 16a: Bulk RNA expression associations

`00_calculate_bulk_tpm_from_counts.R` calculates the deposited TPM table from the featureCounts output and GRCh38 Ensembl release 108 gene-name mapping in `00_Global_Data/Shared_Inputs`. It retains the first entry for each annotated gene name, selects the 16 stimulated PBMC samples and calculates TPM over all 39,453 retained genes. The output is `00_Global_Data/Fetal_Immune_Atlas_PHA_Bulk_TPM.csv`.

`01_test_gene_expression_sample_order_associations.R` reads that TPM table and the sample metadata from `00_Global_Data`. It extracts 61 genes corresponding to Olink assays and applies log2(TPM + 1). The manuscript analysis uses the age-ordered sample positions (1–16) as the Spearman predictor, R's default Spearman test and BH correction across the 61 genes. Expression values and correlation results are saved in `00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_16/a_S16a/03_PlotData`, which is read by both panels.

`02_plot_gene_age_correlation_bars.R` plots the 61 correlation coefficients with the manuscript P-value color gradient. It writes PDF and PNG files to `../02_Figures`.
