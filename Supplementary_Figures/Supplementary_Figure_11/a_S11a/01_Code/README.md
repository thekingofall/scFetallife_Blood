# Supplementary Figure 11a: TCR clonal expansion on the UMAP

`01_extract_tcr_clone_umap_data.py` reads clone annotations and the stored UMAP coordinates from `00_Global_Data/Fetal_Immune_Atlas_TCR.h5ad`. It exports a table to `../03_PlotData`.

`02_plot_tcr_clone_size_umap.R` plots cells in the clone-size groups (1, 2 and at least 3 cells), using gray, blue and red. PDF and PNG files are saved in `../02_Figures`.
