# Figure 2e: Temporal expression programs

Five gene programs are calculated across five post-conception age intervals using stage-average expression profiles from 21 PBMC specimens, with P14 excluded. The intervals contain 2, 5, 5, 4 and 5 specimens, respectively. Clusters 10, 5, 1, 8 and 11 contain 1,422, 876, 1,503, 1,303 and 1,126 genes, respectively. The input files contain the stage-average profiles and Mfuzz module assignments.

`01_define_expression_clustering.R` provides the ClusterGVis 0.1.0 clustering function, including expression standardization and the Mfuzz seed of 5201314.

`02_cluster_pbmc_age_profiles.R` reads `00_Global_Data/Fetal_Immune_Atlas_scRNA_Pseudobulk_TPM.csv`, excludes P14 and calculates the five stage-average profiles. It saves them to `00_Global_Data/Shared_Inputs/Main_Figure_02/e_M02e/M02e_PBMC_Five_Stage_Mean_TPM.csv`, standardizes each gene and fits 16 Mfuzz modules. It writes the assignments and module sizes to `../03_PlotData`.

`03_enrich_temporal_programs.R` maps symbols to Entrez identifiers and runs GO biological-process enrichment for the five selected programs. All tested terms and BH-adjusted P values are saved in `../03_PlotData/PBMC_Program_GO_BP.csv`.

`04_plot_temporal_heatmap.R` reads the module assignments and program labels, orders genes within each program by hierarchical clustering, and draws the heatmap. Gene counts are calculated from the assignments. PDF and PNG files are saved in `../02_Figures`; the ordered expression values are saved in `../03_PlotData`.

Scripts 02 and 03 calculate gene clusters and GO enrichment from the stage-average expression table. Script 04 draws the heatmap from the module assignments.
