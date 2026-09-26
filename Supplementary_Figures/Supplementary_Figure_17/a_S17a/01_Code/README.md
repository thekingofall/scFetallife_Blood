# Supplementary Figure 17: Temporal expression programs in stimulated PBMCs

The heatmap shows six programs from the 16-sample stimulated PBMC bulk RNA-seq cohort. The selected Mfuzz modules are 9, 8, 1, 7, 6 and 5, containing 1,418, 1,704, 1,655, 1,740, 1,734 and 2,556 genes.

`01_cluster_bulk_age_profiles.R` reads gene counts and gene lengths from `00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_17/00_Counts_and_Gene_Lengths/PBMC_Counts_and_Gene_Lengths.rds`. It calculates TPM with the log-space formula, averages each gene within four age groups, and fits 16 Mfuzz modules with seed 5201314 using the shared clustering function in Figure 2e. Age groups follow the whole-week ages used in the manuscript analysis; decimal post-conception ages are rounded to the nearest week for grouping. The sample assignments, stage means and module memberships are saved in `../03_PlotData`.

`00_Global_Data/Fetal_Immune_Atlas_PHA_Bulk_TPM.csv` provides the TPM table for the same 16 samples. The clustering script calculates TPM from the counts and gene lengths.

`02_enrich_bulk_programs.R` runs GO biological-process enrichment for the six selected modules. It saves all tested terms, their gene counts and BH-adjusted P values in `../03_PlotData/PBMC_Bulk_Program_GO_BP.csv`. The heatmap displays the selected GO terms for each module.

`03_plot_bulk_temporal_heatmap.R` orders genes within each module, draws the four-stage heatmap and calculates the displayed gene counts from the module assignments. It saves the ordered plot values in `../03_PlotData` and PDF/PNG figures in `../02_Figures`.

Run scripts 01–03 to repeat the analysis from gene counts, or script 03 to redraw the supplied results.
