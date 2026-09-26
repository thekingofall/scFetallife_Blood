# Supplementary Figure 5c: Cell-type composition across datasets

## Scripts

### `01_compare_fetal_and_cord_blood_cell_composition.py`

Reads the integrated fetal/reference cell metadata and UMAP coordinates. It draws the combined atlas, dataset-specific embeddings, and sample compositions for Supplementary Figure 5a-c. For panel d, it compares late-fetal and cord-blood cell frequencies with Fisher's exact test, plots log2 odds ratios, and exports raw P values and BH-adjusted values.

## Data

**[S5c_annotation_color_mapping.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/c_S05c/03_PlotData/S5c_annotation_color_mapping.csv)**

Harmonized cell-type labels and colors for the sample-composition plot, with fetal and cord-blood label mappings.

**[S5c_composition_complete.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/a_S05a/03_PlotData/S5c_composition_complete.csv)**

Cell counts, sample totals, and cell-type proportions for the fetal/reference comparison, with age, stage, and dataset labels.

**[S5c_composition_complete.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/a_S05a/03_PlotData/S5c_composition_complete.csv)**

Cell counts, sample totals, and cell-type proportions for the fetal/reference comparison, with age, stage, and dataset labels.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
