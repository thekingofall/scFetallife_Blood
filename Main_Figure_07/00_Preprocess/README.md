# Figure 7 preprocessing

`01_prepare_olink_npx.R` reads the `plasma` and `cell supernatant` sheets of `00_Global_Data/Fetal_Immune_Atlas_Olink_NPX.xlsx`. It selects the samples recorded for each assay in `Fetal_Immune_Atlas_Sample_Metadata.csv`, identifies Olink protein columns, matches sample ages and converts the workbook into long tables. NPX values are rounded to two decimal places.

Outputs: `02_Processed/M07_Plasma_NPX.csv` and `02_Processed/M07_Stimulated_NPX.csv`, with `MainID`, `Post_Conception_Age_Weeks`, `Protein` and `NPX` columns.

`Main_Figure_08/00_MOFA_Analysis/01_Code/02_prepare_protein_flow_bulk_views.R` reads the same Olink workbook. Figure 7 plotting scripts read the tables in each panel's `03_PlotData` folder. The preprocessing outputs contain sample-level NPX values and ages for protein analyses.
