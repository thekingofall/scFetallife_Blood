# Multiomics factor analysis

Figure 8 integrates 25 views from 16 fetal blood samples: 18 cell-type expression profiles, stimulated bulk RNA-seq, plasma and stimulated Olink measurements, spectral-flow frequencies, CD4 TCR, CD8 TCR, and BCR features. The 18,168 input features include 80 BCR heavy-chain gene-usage features (IGHV, IGHD and IGHJ). The MOFA model contains 15 active factors.

## Data preparation and analysis

Run scripts 01–04 in order to prepare the model inputs. R scripts use the library configuration in `00_Settings/00_setup_R.R`. Script 03 requires Python with `anndata` and `pandas`.

| Script | Input and processing | Output in `02_Data` |
|---|---|---|
| `01_prepare_single_cell_views.R` | Selects blood samples with all five assay availability flags in the project metadata. Reads expression sums grouped by sample and cell type, together with the selected gene list. Equalizes library totals within each cell type, retains the selected genes, applies `log2(x + 1)` and quantile normalization, and fills absent sample–cell-type combinations with zero. | `MOFA_Sample_Metadata.csv`, `Single_Cell_Features.csv.gz` |
| `02_prepare_protein_flow_bulk_views.R` | Reads both sheets of the project Olink workbook and retains protein assays with Olink IDs. Stores NPX to two decimals using the source-table rounding convention. Reads the complete parent-gate flow-frequency table. Filters stimulated bulk TPM by a summed TPM above 10, applies `log2(TPM + 1)`, and retains genes whose transformed values sum to more than 100 across samples. | `Protein_Flow_Bulk_Features.csv.gz` |
| `03_prepare_repertoire_views.py` | Selects the 80 IGH gene-usage features for the BCR view. Combines CD4 and CD8 TCR gene usage with CDR3 length frequencies calculated from the project TCR h5ad metadata. Length counts are divided by the total across the selected samples within each T-cell view. | `Repertoire_Features.csv.gz` |
| `04_assemble_mofa_views.R` | Joins the three feature tables by sample, retains missing values, and applies an inverse-normal rank transformation to each feature across samples. Removes features matching `__RPS`, `__RPL` or `^RP`, separates the single-cell profiles into cell-type views, and aligns every matrix to the same sample order. | `MOFA_Normalized_Features.csv.gz`, `MOFA_View_Inputs.rds`, `MOFA_View_Summary.csv` |
| `05_fit_multiomics_factors.R` | Fits MOFA to the matrices produced by script 04, with 16 initial factors, seed 2024, view scaling, and up to 50,000 iterations. Requires MOFA2 and its `mofapy2` Python dependency. | `MOFA_Model.hdf5` |
| `06_extract_factor_results.R` | Extracts explained variance, sample factor scores, factor–age Spearman correlations with BH correction, and positive and negative Factor 3 feature weights. Exports the tables for Figure 8b–f and Supplementary Figure 18. | Factor and feature tables; panel-specific tables under `03_PlotData` |
| `07_enrich_factor3_features.R` | Runs GO biological-process enrichment for the 100 strongest positive and negative Factor 3 features in each of the 19 gene-expression views. Uses SYMBOL identifiers, BH correction, minimum term size 1, and P-value and q-value cutoffs of 0.05. | `MOFA_Factor3_Top100_Features.csv`, `MOFA_Factor3_GO_BP.csv`, `MOFA_Factor3_Selected_GO_BP.csv` |
| `08_assemble_figure8.R` | Combines the manuscript integration schematic and the five computed panels into Figure 8. Requires `png`, `grid`, `patchwork` and `ragg`. | Complete PDF and PNG in `02_Figures` |

## Inputs

Primary files are read from `00_Global_Data`: the sample metadata, Olink workbook, stimulated bulk TPM, and TCR h5ad file. These are the project copies of the deposited data.

`00_Global_Data/Shared_Inputs/Main_Figure_08/00_MOFA/01_Source_Data` contains cell-type expression sums, the 60-feature flow table, receptor-gene usage and CDR3 length definitions. `Single_Cell_Selected_Genes.csv` specifies the genes included in each cell-type view.

## Figures

To draw Figure 8 and Supplementary Figure 18, run script 06 and then the plotting scripts in `b_M08b` through `f_M08f` and `Supplementary_Figures/Supplementary_Figure_18/a_S18a`. Script 06 reads `00_Global_Data/Shared_Inputs/Main_Figure_08/00_MOFA/Fetal_Immune_Atlas_MOFA.hdf5`.

Script 05 saves its fitted model to `02_Data/MOFA_Model.hdf5`. Scripts 06 and 07 read `00_Global_Data/Shared_Inputs/Main_Figure_08/00_MOFA/Fetal_Immune_Atlas_MOFA.hdf5`. The GO plotting tables are exported by script 06 from `MOFA_GO_BP_All.csv` and `MOFA_GO_BP_Selected.csv`; script 07 calculates GO enrichment from the model feature weights.
