# Figure 2d: Hotspot module enrichment

`01_assign_hotspot_modules_and_run_go.R` assigns the early and late Hotspot gene modules to the reference lineages by gene overlap and runs GO enrichment with clusterProfiler. GO-BP, GO-CC and GO-MF are tested with BH correction and the default annotated-gene background. Results are saved under `02_Data/02_Processed`.

`02_prepare_selected_go_terms.R` reads the 24 terms selected for the manuscript from `03_PlotData/M02d_Selected_GO_Terms.xlsx`. It selects the P values, supporting genes and gene counts, orders the terms within each module and stage, and writes `03_PlotData/M02d_selected_pathways.csv`.

`03_plot_selected_go_terms.R` draws the six module panels. Dot size shows the supporting gene count; color shows -log10(BH-adjusted P), capped at 12. Each module contains two early and two late terms. The two early NK terms are GO-CC; the other 22 terms are GO-BP. PNG and PDF files are written to `02_Figures`.

Scripts 01 and 02 calculate enrichment and prepare the selected terms. Script 03 draws the panel from `03_PlotData/M02d_selected_pathways.csv`.
