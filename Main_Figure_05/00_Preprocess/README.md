# Figure 5 preprocessing

| Script | Processing | Output |
| --- | --- | --- |
| `01_merge_bcr_contigs_and_filter_pairs.py` | Formats productive-contig annotations, matches sample-prefixed barcodes to the atlas and assigns chain-pairing classes. Retains B cells with a single pair or an extra VJ/VDJ chain and applies the curated BCR cell list. | `02_Processed/M05_BCR_Contigs.csv`, `../../00_Global_Data/Fetal_Immune_Atlas_BCR.h5ad` |
| `02_assign_bcr_clonotypes.py` | Converts the BCR object to Dandelion, retains annotated V/J and CDR3 sequences, defines clones with the `hh_s5f` distance threshold 0.15858785672187678, and calculates the clone network and capped clone sizes. | `02_Processed/M05_BCR_Clone_Annotations.csv`, `M05_BCR_Clonotypes.pkl`, `M05_BCR_Clonotypes.h5ad` |

Input contigs are in `00_Global_Data/Shared_Inputs/Main_Figure_04_05/00_Receptor_Annotations/01_Aggregated_Contigs/Gao_All_BCR.csv`. `01_Inputs/BCR_Cells.csv.gz` contains the 30,173 cells in the deposited BCR object. The paired-receptor h5ad is written by step 1; step 2 writes the Dandelion clone-analysis outputs.

Scirpy 0.12.0 processes the `IR_*` receptor annotations. Dandelion clone calling also requires its R/Shazam and graph-layout dependencies.
