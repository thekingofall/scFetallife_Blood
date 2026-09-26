# Figure 6 preprocessing

`01_prepare_flow_events.py` runs the flow preprocessing in `../a_M06a/01_Code/01_gate_flow_events_and_compute_umap.py`.

It reads the FlowJo workspaces and FCS events, applies the gates, selects the recorded downsampled populations, compensates and transforms marker intensities, and calculates two-dimensional UMAP coordinates with 15 neighbors and `min_dist=0.5`.

Download the five FCS ZIP files and three workspace/metadata files from [Zenodo](https://doi.org/10.5281/zenodo.22952935). Follow the [flow-data setup instructions](../../00_Settings/README.md) to place them under `00_Global_Data/Flow_Data/fetal flow data`.

The Figure 6 plotting code reads `M06a_S14_PBMC_event_coordinates.csv.gz` and `M06a_S14_Organ_event_coordinates.csv.gz` from `00_Global_Data/Shared_Inputs/Main_Figure_06/a_M06a`. Sample metadata and transformed marker matrices are written to `00_Global_Data/Shared_Inputs/Main_Figure_06/a_M06a`. Figure 6 uses flow-cytometry data and has no h5ad input.
