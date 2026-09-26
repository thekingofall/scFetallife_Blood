# Figure 5c: BCR chain pairing

`01_classify_receptor_chain_pairing.py` reads the cell atlas and combined BCR contig annotations from `00_Global_Data`. It joins contigs to atlas cells by barcode and classifies receptor chains with Scirpy `chain_qc`. Cells with a receptor annotation are retained with their atlas UMAP coordinates. Aligned contigs and cell-level plotting data are written to `../03_PlotData`.

`02_plot_receptor_chain_pairing_umap.R` draws the chain-pairing classes using the manuscript colors and exports PDF and PNG files to `../02_Figures`. Colors use the same categories as panels a and b.

Run the Python script first, then the R script. The chain-classification analysis uses Scirpy 0.12.0 and AnnData 0.10.7. Input paths are resolved from the script location.
