# Figure 8 preprocessing

| Script | Processing | Output |
| --- | --- | --- |
| `01_aggregate_single_cell_expression.py` | Selects PBMC samples with all five assays and the 18 MOFA cell-type views. Sums the atlas log-normalized expression within each sample/cell-type group. | `00_Global_Data/Shared_Inputs/Main_Figure_08/00_MOFA/01_Source_Data/Single_Cell_Pseudobulk.csv.gz` |
| `02_prepare_multiomics_views.py` | Runs the RNA, protein/flow/bulk, receptor-feature and view-assembly scripts in order. These normalize the features, calculate TCR CDR3-length frequencies, select BCR gene usage, align the samples and assemble the MOFA views. | `../00_MOFA_Analysis/02_Data/Single_Cell_Features.csv.gz`, `Protein_Flow_Bulk_Features.csv.gz`, `Repertoire_Features.csv.gz`, `MOFA_Normalized_Features.csv.gz`, `MOFA_View_Inputs.rds` |

The first script aggregates log-normalized expression by sample and cell type. The second runs the preparation scripts in `../00_MOFA_Analysis/01_Code` using the project R configuration. Set `RSCRIPT` when Rscript is outside the executable search path.

The h5ad inputs are produced in the Figure 1 and Figure 4 preprocessing folders. The fitted model used for Figure 8 is `00_Global_Data/Shared_Inputs/Main_Figure_08/00_MOFA/Fetal_Immune_Atlas_MOFA.hdf5`. Model fitting and factor-result extraction scripts are in `../00_MOFA_Analysis/01_Code`.
