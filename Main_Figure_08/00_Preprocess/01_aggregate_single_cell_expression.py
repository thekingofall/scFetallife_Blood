"""Sum atlas log-normalized expression by sample and cell type for the MOFA input."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "00_Settings"))

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from preprocess_inputs import GLOBAL
atlas = ad.read_h5ad(GLOBAL / "Fetal_Immune_Atlas.h5ad", backed="r")
metadata = pd.read_csv(GLOBAL / "Fetal_Immune_Atlas_Sample_Metadata.csv")
assays = ["scRNAseq_Available", "Flow_Cytometry_Available", "Plasma_Olink_Available", "Stimulated_Olink_Available", "PHA_Bulk_RNA_Available"]
samples = metadata.loc[metadata["Main_Organ"].eq("PBMC") & metadata[assays].eq(1).all(axis=1), "MainID"]
selection = pd.read_csv(GLOBAL / "Shared_Inputs/Main_Figure_08/00_MOFA/01_Source_Data/Single_Cell_Selected_Genes.csv")
obs = atlas.obs.loc[atlas.obs["MainID"].astype(str).isin(samples) & atlas.obs["Last_cell_type"].astype(str).isin(selection["cell_type"])].copy()
expression = atlas[obs.index].to_memory()
keys = obs["MainID"].astype(str) + "_" + obs["Last_cell_type"].astype(str)
codes, groups = pd.factorize(keys, sort=True)
summing = sparse.csr_matrix((np.ones(len(obs)), (codes, np.arange(len(obs)))), shape=(len(groups), len(obs)))
values = summing @ expression.X
output = GLOBAL / "Shared_Inputs/Main_Figure_08/00_MOFA/01_Source_Data/Single_Cell_Pseudobulk.csv.gz"
pd.DataFrame(values.toarray(), index=groups, columns=expression.var_names).to_csv(output, index_label="sample_cell_type")
atlas.file.close()
