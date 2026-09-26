# Supplementary Figure 5d: Late-fetal versus cord-blood cell enrichment

## Scripts

### `01_compare_fetal_and_cord_blood_cell_composition.py`

Reads the integrated fetal/reference cell metadata and UMAP coordinates. It draws the combined atlas, dataset-specific embeddings, and sample compositions for Supplementary Figure 5a-c. For panel d, it compares late-fetal and cord-blood cell frequencies with Fisher's exact test, plots log2 odds ratios, and exports raw P values and BH-adjusted values.

## Data

**[S5d_contingency.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/a_S05a/03_PlotData/S5d_contingency.csv)**

Late-fetal and cord-blood cell counts for each cell type, used to construct the Fisher-test contingency tables.

**[S5d_fisher_rawP_BHsecondary.csv](../../../../00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/a_S05a/03_PlotData/S5d_fisher_rawP_BHsecondary.csv)**

Late-fetal versus cord-blood enrichment results: cell counts, odds ratios, log2 odds ratios with a 0.5 pseudocount, raw P values, and BH-adjusted values.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
