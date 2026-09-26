# Figure 1b: sample overview data

`M01b_multimodal_plot_data.csv` contains the 37 sample and tissue entries plotted in Figure 1b. `MainID` matches the shared sample metadata. `AlmostWeek` records the manuscript display positions, including offsets separating nearby samples. The five assay flags describe the assays included in the manuscript and are refreshed from `00_Global_Data/Fetal_Immune_Atlas_Sample_Metadata.csv` when the plotting script runs.

`M01b_cell_counts.tsv` contains the scRNA-seq cell total for each tissue.

`M01b_sample_overview_plot_data.tsv` contains single-cell sample ages and scRNA-seq, TCR and BCR cell counts.

`data_dependencies.tsv` lists the input tables and plotting script.
