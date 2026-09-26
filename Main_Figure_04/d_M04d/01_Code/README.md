# Figure 4d: UMAP of TCR cells by tissue

## Scripts

### `01_plot_receptor_umap_by_tissue.py`

Reads the existing TCR and BCR UMAP coordinates and colors cells by `Main_Organ`. It draws Figures 4d and 5d using the four tissue colors and writes PDF/PNG figures plus the cell-level plotting tables. Set `CLEAN_ROOT` to choose the output root.

## Data

**[M04c_TCR_Chain_Pairing_UMAP.csv](../../../00_Global_Data/Shared_Inputs/Main_Figure_04_05/01_Chain_Pairing_Clonotypes/M04c_TCR_Chain_Pairing_UMAP.csv)**

Existing UMAP coordinates for receptor-annotated cells, with chain-pairing status, cell types, samples, and tissues.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../00_Settings/00_External_Data_Dependencies.tsv).
