# Supplementary Figure 5b: UMAP by dataset

## Scripts

### `01_compare_fetal_and_cord_blood_cell_composition.py`

Reads the integrated fetal/reference cell metadata and UMAP coordinates. It draws the combined atlas, dataset-specific embeddings, and sample compositions for Supplementary Figure 5a-c. For panel d, it compares late-fetal and cord-blood cell frequencies with Fisher's exact test, plots log2 odds ratios, and exports raw P values and BH-adjusted values.

## Data

**[S5b_source_mapping.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/b_S05b/03_PlotData/S5b_source_mapping.csv)**

Dataset names, display order, and foreground/background colors for the dataset-specific UMAP panels.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
