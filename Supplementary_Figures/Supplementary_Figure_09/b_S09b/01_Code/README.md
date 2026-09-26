# Supplementary Figure 9b: Incoming signaling to NK T cells

`01_extract_incoming_interactions.R` reads the saved early- and late-stage PBMC communication networks from `00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_09/01_CellChat_Models/S09_CellChat_Communication_Networks.rds`. It extracts interactions targeting NK T cells and selects the pathways shown in the manuscript. The networks contain fitted interaction probabilities, permutation P values and cell-type annotations.

The plotted communication probability follows CellChat's `-1/log(probability)` transformation. Significant interactions from either stage determine the displayed ligand-receptor rows. The table and its factor order are saved in `../03_PlotData`.

`02_plot_incoming_interactions.R` draws early and late stages in separate panels. Rows show ligand-receptor pairs, columns show sender cell types, color shows transformed probability and point size shows permutation significance. It writes PDF and PNG files to `../02_Figures`.

Run scripts 01 and 02 in order to redraw the panel from the supplied communication networks.
