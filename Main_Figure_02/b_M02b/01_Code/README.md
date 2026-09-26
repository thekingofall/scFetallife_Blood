# Figure 2b: Age-associated changes in PBMC composition

## Scripts

### `01_analyze_cell_proportion_age_associations.R`

Reads the PBMC cell-composition and erythroid summary tables and relates cell proportions to `Post_Conception_Age_Weeks`. It calculates Spearman correlations, draws the selected age-association panels, and exports the plotted values, full statistics, significant results, and PDF/PNG figures.

## Data

**`M02b_PBMC_age_correlations_plot_data.csv`**

Sample-level PBMC cell-type counts and proportions, with post-conception age and the display labels used for the age-association plots.

**[M02b_PBMC_age_correlations_significant.csv](../../../00_Global_Data/Shared_Inputs/Main_Figure_02/b_M02b/03_PlotData/M02b_PBMC_age_correlations_significant.csv)**

Spearman cell-proportion/age results, including sample counts, raw P values, q values, and significance flags. This file contains the selected significant results.

**[M02b_PBMC_age_correlations_statistics.csv](../../../00_Global_Data/Shared_Inputs/Main_Figure_02/b_M02b/03_PlotData/M02b_PBMC_age_correlations_statistics.csv)**

Spearman cell-proportion/age results, including sample counts, raw P values, q values, and significance flags.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../00_Settings/00_External_Data_Dependencies.tsv).
