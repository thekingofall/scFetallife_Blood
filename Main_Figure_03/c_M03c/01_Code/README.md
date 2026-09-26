# Figure 3c: Changes in incoming signaling

`01_plot_incoming_signaling_changes.R` reads [M03c_signaling_changes_data.csv](../../../00_Global_Data/Shared_Inputs/Main_Figure_03/c_M03c/03_PlotData/M03c_signaling_changes_data.csv) and plots pathway-level changes for CX3CR1+ NK, CXCR6+ NK, CD56highCD16low NK, naive CD8 T and NK T cells. The shared CellChat analysis retains incoming changes above 0.3 for the first four populations and below -0.3 for NK T cells.

Run `Main_Figure_03/00_CellChat_Analysis/01_Code/02_summarize_interaction_strength_and_signaling.R` to recalculate this table from the early- and late-stage models. PDF and PNG files are saved in `../02_Figures`.
