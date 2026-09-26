# Supplementary Figure 5a: Integrated fetal and reference blood atlas

## Scripts

### `01_compare_fetal_and_cord_blood_cell_composition.py`

Reads the integrated fetal/reference cell metadata and UMAP coordinates. It draws the combined atlas, dataset-specific embeddings, and sample compositions for Supplementary Figure 5a-c. For panel d, it compares late-fetal and cord-blood cell frequencies with Fisher's exact test, plots log2 odds ratios, and exports raw P values and BH-adjusted values.

## Data

**[S5a_annotation_color_mapping.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/a_S05a/03_PlotData/S5a_annotation_color_mapping.csv)**

Integrated-atlas cell labels and colors, including the correspondence to cord-blood labels.

**[S5ab_cell_metadata_coordinates.csv.gz](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/a_S05a/03_PlotData/S5ab_cell_metadata_coordinates.csv.gz)**

Integrated fetal/reference UMAP coordinates, dataset and stage labels, and the cell-type annotations used by the comparison panels.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
