# Figure 2c: Early and late gene-correlation modules

`01_plot_gene_correlation_heatmaps.R` reads the early- and late-stage Hotspot correlation matrices and module assignments from `00_Global_Data/Shared_Inputs/Main_Figure_02/c_M02c`. The table `../03_Source_Data/M02c_fresh_to_canonical_module_mapping.csv` assigns numerical modules to seven lineage programs.

Genes are ordered within each program by hierarchical clustering of the local-correlation matrix. The script draws the early and late heatmaps and exports the gene order and program sizes. PDF and PNG files are saved in `../02_Figures`; module counts are saved in `../03_PlotData`.

GO enrichment for these programs is provided in the Figure 2d folder.
