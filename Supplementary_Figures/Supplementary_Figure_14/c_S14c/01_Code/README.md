# Supplementary Figure 14c: PBMC flow-cytometry UMAP by sample

## Scripts

### `01_gate_flow_events_and_compute_umap.py`

Reads the blood and organ FCS files, the FlowJo workspaces, and sample metadata. It replays compensation, transformations, and gates with FlowKit, applies the downsampling rules specified in the FlowJo workspaces, and computes UMAP coordinates. It saves event coordinates in `00_Global_Data/Shared_Inputs/Main_Figure_06/a_M06a`. Sample metadata and transformed marker measurements are saved in `00_Global_Data/Shared_Inputs/Main_Figure_06/a_M06a`. Each run calculates the gated events from the FCS files. Events matching several terminal gates take the last matching category in the cell-type order.

### `02_plot_flow_umap.py`

Reads the blood and organ event-coordinate tables in `00_Global_Data/Shared_Inputs/Main_Figure_06/a_M06a` and draws the combined flow-cytometry UMAP, tissue panels, and individual-sample panels. It controls point sizes, cell-type labels, legends, and panel layout, and exports PDF/PNG figures. Gate assignments and coordinates come from the input tables.

## Data

**[M06a_S14_PBMC_event_coordinates.csv.gz](../../../../00_Global_Data/Shared_Inputs/Main_Figure_06/a_M06a/M06a_S14_PBMC_event_coordinates.csv.gz)**

Selected flow-cytometry events with sample, tissue, cell-type gate, downsampling information, and the two UMAP coordinates.

**[M06a_S14_PBMC_transformed_marker_values.csv.gz](../../../../00_Global_Data/Shared_Inputs/Main_Figure_06/a_M06a/M06a_S14_PBMC_transformed_marker_values.csv.gz)**

Compensated and transformed flow-cytometry marker measurements for each selected event. Join to the coordinate table by `sample_id` and `event_index`.

**[M06a_S14_sample_inventory.csv](../../../../00_Global_Data/Shared_Inputs/Main_Figure_06/a_M06a/M06a_S14_sample_inventory.csv)**

Flow sample identifiers, parent-population counts, downsampling targets, selected-event counts, and marker-panel information.

Shared and upstream input paths are listed in `../03_PlotData/data_dependencies.tsv` and the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
