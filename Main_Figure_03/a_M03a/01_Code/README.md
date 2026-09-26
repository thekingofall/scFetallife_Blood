# Figure 3a: Incoming and outgoing cell-contact signaling

## Scripts

### `01_plot_cellchat_signaling_roles.R`

Reads the signaling-role tables and draws the five tissue/stage panels in Figure 3a. Outgoing and incoming interaction strength determine the point coordinates; network link counts determine point size. It exports the plotting values and PDF/PNG figures.

## Running the analysis

Run `01_plot_cellchat_signaling_roles.R` to read the shared signaling-role tables and draw the five tissue/stage panels.

## Data

**`Figure3a_signaling_roles_plot_source.csv`**

Incoming and outgoing signaling strengths for each cell type and tissue/stage group, with network link counts and label/point-size settings.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../00_Settings/00_External_Data_Dependencies.tsv).
