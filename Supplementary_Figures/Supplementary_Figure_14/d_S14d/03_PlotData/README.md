# Supplementary Figure 14d: Organ flow-cytometry UMAP by sample

**[M06a_S14_Organ_event_coordinates.csv.gz](../../../../00_Global_Data/Shared_Inputs/Main_Figure_06/a_M06a/M06a_S14_Organ_event_coordinates.csv.gz)**

Selected flow-cytometry events with sample, tissue, cell-type gate, downsampling information, and the two UMAP coordinates.

**[M06a_S14_Organ_transformed_marker_values.csv.gz](../../../../00_Global_Data/Shared_Inputs/Main_Figure_06/a_M06a/M06a_S14_Organ_transformed_marker_values.csv.gz)**

Compensated and transformed flow-cytometry marker measurements for each selected event. Join to the coordinate table by `sample_id` and `event_index`.

**[M06a_S14_sample_inventory.csv](../../../../00_Global_Data/Shared_Inputs/Main_Figure_06/a_M06a/M06a_S14_sample_inventory.csv)**

Flow sample identifiers, parent-population counts, downsampling targets, selected-event counts, and marker-panel information.

See `../01_Code/README.md` for the scripts that use or generate these tables. `data_dependencies.tsv` records shared-resource paths; additional inputs are listed in the project-level [input-file list](../../../../00_Settings/00_External_Data_Dependencies.tsv).
