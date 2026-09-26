# Expression datasets

`Fetal_Immune_Atlas_scRNA_Pseudobulk_TPM.csv` contains gene-level expression aggregated from single-cell RNA-seq for 21 PBMC samples (20,719 genes). Figure 2e reads the five-stage expression means in `Shared_Inputs/Main_Figure_02/e_M02e`.

`Fetal_Immune_Atlas_PHA_Bulk_TPM.csv` contains TPM values from stimulated PBMC bulk RNA-seq for 16 samples (39,453 genes). It supplies Supplementary Figures 16 and 17. Columns are MainID sample identifiers. Values are TPM, before log transformation; each sample sums to 1,000,000. Gene annotation uses GRCh38 Ensembl release 108; duplicated gene names retain the first entry.

`Fetal_Immune_Atlas_Sample_Metadata.csv` records availability of the two assays in separate columns. The bulk TPM calculation is provided in `Supplementary_Figures/Supplementary_Figure_16/a_S16a/01_Code/00_calculate_bulk_tpm_from_counts.R`; the figure analysis reads the supplied TPM table.
