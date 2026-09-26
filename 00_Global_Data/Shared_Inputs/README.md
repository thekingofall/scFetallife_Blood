# Shared analysis inputs

The atlas, receptor data and deposited assay tables are in the parent `00_Global_Data` directory. These folders contain annotations, fitted models and intermediate data used by the manuscript figures.

| Folder | Contents |
| --- | --- |
| `00_Common` | Sample metadata, GRCh38 gene annotations and receptor gene order. |
| `Main_Figure_01` | Atlas cell annotations, UMAP coordinates and cell-type colors; also used by Figures 2 and S1-S3. |
| `Main_Figure_02` | Cell proportions, Hotspot modules, gene correlations, GO results and five-stage expression averages. |
| `Main_Figure_03` | CellChat models, signaling roles and ligand-receptor summaries. |
| `Main_Figure_04_05` | TCR/BCR contig annotations, chain pairing, gene usage, VJ/DJ pairing and CDR3 lengths. Several scripts process both receptors together. |
| `Main_Figure_06` | Flow-cytometry summaries and event coordinates; also used by Figures S14 and S15. |
| `Main_Figure_07` | Plasma and stimulation protein measurements, age-association statistics and xCell scores. |
| `Main_Figure_08` | MOFA source data, 25-view input matrices, the fitted model and GO enrichment results. |
| `Supplementary_Figures` | Inputs grouped by supplementary figure number. |

Within each figure folder, panel labels follow the manuscript order. `M` denotes a main figure and `S` a supplementary figure. Sample-level contig filenames retain their sample identifiers.

Figure 2e reads five-stage expression averages and gene-program assignments based on 22 PBMC specimens. These tables are in `Main_Figure_02/e_M02e`. The deposited single-cell pseudobulk table contains 21 specimens; this panel uses the stage-average input tables.

Figure S9 uses the early- and late-stage communication networks under `Supplementary_Figures/Supplementary_Figure_09/01_CellChat_Models`. Figure 3 uses the 28-cell-type models [PBMC_Early_exact28_cellchat.rds](Main_Figure_03/00_CellChat_Models/PBMC_Early_exact28_cellchat.rds) and [PBMC_Late_exact28_cellchat.rds](Main_Figure_03/00_CellChat_Models/PBMC_Late_exact28_cellchat.rds) in `Main_Figure_03/00_CellChat_Models`.

Figure 8 uses 16 samples and 25 views; the manuscript model has 15 active factors. The preprocessing, fitting and enrichment scripts are in `Main_Figure_08/00_MOFA_Analysis` at the repository root. See [MOFA source data](Main_Figure_08/00_MOFA/01_Source_Data/README.md) for its inputs.

Figure S17 uses the counts and gene lengths in `Supplementary_Figures/Supplementary_Figure_17/00_Counts_and_Gene_Lengths`. The analysis calculates TPM and fits temporal modules from four stage averages.

The repository input file list is in `00_Settings/00_External_Data_Dependencies.tsv`.

Plotting tables shared by multiple scripts are stored once in their figure and panel subdirectories. Panel READMEs link to the tables they use.
