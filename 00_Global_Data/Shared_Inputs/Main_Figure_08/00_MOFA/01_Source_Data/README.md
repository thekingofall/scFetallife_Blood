# MOFA source data

These tables supply the preprocessing steps in `Main_Figure_08/00_MOFA_Analysis/01_Code`.

| File | Contents |
|---|---|
| `Single_Cell_Pseudobulk.csv.gz` | Gene-expression sums for observed sample–cell-type combinations in the 18 single-cell views. Row identifiers contain the sample and cell type. Values are sums of atlas log-normalized expression within each sample and cell type. |
| `Single_Cell_Selected_Genes.csv` | Cell-type-specific genes included in the MOFA analysis. Columns are `cell_type` and `gene`. |
| `Flow_Parent_Frequencies.csv` | Sixty flow-cytometry measurements expressed as percentages of the corresponding parent gates, with one sample per column. This complete integration table includes gates beyond the cell populations listed in the deposited proportions workbook. |
| `Receptor_Gene_Usage.csv` | Sample-level receptor-gene usage for the CD4 TCR, CD8 TCR, and BCR views. Columns are `MainID`, `view`, `feature`, and `value`. Values are receptor-gene usage percentages. Figure 8 selects IGHV, IGHD and IGHJ genes for its 80-feature BCR view. |
| `Receptor_CDR3_Features.csv` | CDR3 length-bin definitions. Figure 8 uses the CD4 TCR and CD8 TCR bins, whose frequencies are computed from the project TCR h5ad metadata. Length bins without observations in a sample receive zero. |

The integration sample set is selected from `Fetal_Immune_Atlas_Sample_Metadata.csv`. Olink NPX and stimulated bulk TPM are read from the primary files in `00_Global_Data`.
