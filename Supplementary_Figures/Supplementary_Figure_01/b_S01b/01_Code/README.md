# Supplementary Figure 1b: Weighted cell abundance across tissues

## Scripts

### `01_plot_age_weighted_cell_abundance.R`

Counts cells by sample and cell type from the shared annotation table, calculates the age-group-weighted abundance for each organ, and draws the Supplementary Figure 1b heatmap. It exports per-sample counts, weighted abundance values, and PDF/PNG figures.

## Data

**[S01b_per_sample_counts_source.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_01/b_S01b/03_PlotData/S01b_per_sample_counts_source.csv)**

Cell-type counts, sample totals, age groups, and proportions used to calculate weighted abundance.

**[S01b_weighted_cell_abundance_source.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_01/b_S01b/03_PlotData/S01b_weighted_cell_abundance_source.csv)**

Weighted cell abundance for each tissue and cell type, with the corresponding lineage.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
