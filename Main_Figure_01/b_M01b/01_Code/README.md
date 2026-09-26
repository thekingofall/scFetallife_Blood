# Figure 1b: samples and assay coverage

`01_plot_sample_overview_and_cell_counts.R` draws the sample overview for blood, liver, thymus and spleen. Assay inclusion is read from `00_Global_Data/Fetal_Immune_Atlas_Sample_Metadata.csv`, the same metadata file supplied with the Zenodo dataset. Samples are matched by `MainID`.

`../03_PlotData/M01b_multimodal_plot_data.csv` supplies the displayed sample positions (`AlmostWeek`). The script updates its assay flags from the shared sample metadata before plotting. `M01b_cell_counts.tsv` supplies the scRNA-seq cell totals for each tissue.

The manuscript includes stimulated proteomics and bulk RNA-seq for 16 PBMC samples. Organ samples are represented by their scRNA-seq and spectral flow cytometry assays. The four tissue rows contain 235,494 cells, including 162,092 blood cells.

Run with `Rscript 01_plot_sample_overview_and_cell_counts.R`. The script loads the shared R settings and saves PDF and PNG files in `../02_Figures`.
