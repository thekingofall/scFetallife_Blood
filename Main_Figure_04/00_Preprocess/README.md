# Figure 4 preprocessing

| Script | Processing | Output |
| --- | --- | --- |
| `01_merge_tcr_contigs_and_filter_pairs.py` | Prefixes contig barcodes with their sample identifiers, merges TCR annotations with the atlas, assigns receptor and chain-pairing classes, removes unpaired and multichain cells, and applies the curated TCR cell list. | `02_Processed/M04_TCR_Contigs.csv`, `M04_Paired_TCR.h5ad` |
| `02_assign_tcr_clonotypes.py` | Defines exact nucleotide clonotypes using both receptor arms and primary chains. Calculates amino-acid alignment clusters with cutoff 15, including a same-V-gene grouping, clone sizes and the clonotype network. | `../../00_Global_Data/Fetal_Immune_Atlas_TCR.h5ad` |

Run the scripts in order. Input contigs are in `00_Global_Data/Shared_Inputs/Main_Figure_04_05/00_Receptor_Annotations/01_Aggregated_Contigs/Gao_All_TCR.csv`. `01_Inputs/TCR_Cells.csv.gz` contains the 53,955 cells in the deposited TCR object. Atlas expression and UMAP coordinates come from `Fetal_Immune_Atlas.h5ad`.

The receptor preprocessing API is Scirpy 0.12.0; specified versions are listed in `../../00_Settings/Preprocess_Python_Versions.tsv`. Clone identifiers and network fits are generated when step 2 runs.
