# Figure 5b: BCR chain pairing by tissue

## Scripts

### `01_define_repertoire_palettes.R`

Defines the color vectors used for tissues, cell populations, receptor-pairing categories, and repertoire plots. Source this file before the R plotting scripts that use these palettes.

### `02_plot_receptor_pairing_by_sample_and_tissue.R`

Reads the TCR and BCR cell-metadata/UMAP tables and tabulates chain-pairing categories by sample and organ. It draws the four stacked-bar panels in Figures 4a-b and 5a-b and writes their proportions and PDF/PNG figures.

## Data

**[M05b_BCR_pairing_by_organ_data.csv](../../../00_Global_Data/Shared_Inputs/Main_Figure_05/b_M05b/03_PlotData/M05b_BCR_pairing_by_organ_data.csv)**

TCR/BCR chain-pairing proportions for each tissue.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../00_Settings/00_External_Data_Dependencies.tsv).
