# Supplementary Figure S11c: TCR clonotype network

`01_plot_tcr_clone_network_by_organ.py` reads `00_Global_Data/Fetal_Immune_Atlas_TCR.h5ad`, orders the stored clonotype-index mapping numerically, and calculates the exact-clonotype network with `min_cells=2` and `random_state=42`. Nodes are colored by tissue using the figure palette. Network coordinates are exported to `../03_PlotData`; figures are saved as PDF and PNG in `../02_Figures`.

The analysis uses Scirpy 0.12.0, Scanpy 1.9.3, AnnData 0.10.7, NumPy 1.23.4 and igraph 0.10.4. Scirpy supplies the network layout and plotting method.
