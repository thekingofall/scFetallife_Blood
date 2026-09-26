# Figure 8: Multiomics integration of fetal blood development

The analysis integrates 25 data views from 16 samples, comprising 18,168 features. The BCR view contains 80 IGHV, IGHD and IGHJ gene-usage features. The supplied MOFA model has 15 active factors. Factor 3 is associated with post-conception age (Spearman R = 0.77353, P = 0.000683, BH-adjusted q = 0.010245).

[Data preparation and model analysis](00_MOFA_Analysis/README.md) describes scripts 01–07. Script `08_assemble_figure8.R` generates the complete [PDF](00_MOFA_Analysis/02_Figures/Figure_08.pdf) and [PNG](00_MOFA_Analysis/02_Figures/Figure_08.png).

| Panel | Contents |
|---|---|
| [a](a_M08a) | Manuscript diagram of the integration framework |
| [b](b_M08b) | Explained variance and feature counts by view |
| [c](c_M08c) | Factor–age correlations with BH correction |
| [d](d_M08d) | Factor 3 scores across post-conception age |
| [e](e_M08e) | The two strongest positive and negative Factor 3 features in each view |
| [f](f_M08f) | GO enrichment of positive and negative Factor 3 features |

Each computed panel reads its compact tables from `03_PlotData` and saves PDF and PNG files in `02_Figures`. The fitted model, normalized input matrices and manuscript GO results are under `../00_Global_Data/Shared_Inputs/Main_Figure_08/00_MOFA`.

## Preprocessing

See [preprocessing scripts and inputs](00_Preprocess/README.md) for the processing order and output filenames.
