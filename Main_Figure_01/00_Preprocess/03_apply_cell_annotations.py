"""Combine normalized expression with the manuscript cell annotations and UMAP coordinates."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "00_Settings"))

import numpy as np
import pandas as pd
import scanpy as sc
from preprocess_inputs import GLOBAL, INPUTS, save_h5ad

processed = Path(__file__).resolve().parent / "02_Processed"
data = sc.read_h5ad(processed / "M01_Filtered_Clustered.h5ad")
annotations = pd.read_csv(GLOBAL / "Shared_Inputs/Main_Figure_01/01_Cell_Annotation/M01_Cell_Annotations_and_UMAP.csv", index_col="Cellname", low_memory=False)
missing = annotations.index.difference(data.obs_names)
if len(missing):
    raise ValueError(f"{len(missing)} curated cells are absent after filtering; check the raw matrices and QC environment")
data = data[annotations.index].copy()
normalized = sc.read_h5ad(processed / "M01_CellTypist.h5ad")
result = normalized[annotations.index].copy()
result.obs = data.obs.copy()
result.obsm = data.obsm.copy()
result.obsp = data.obsp.copy()
result.uns = data.uns.copy()
result.varm.clear()
genes = pd.read_csv(INPUTS / "Atlas_Genes.csv.gz", index_col=0)
result = result[:, genes.index].copy()
result.var = genes
columns = ["Name", "MainID", "Main_Organ", "Gestational_Age_Weeks", "Post_Conception_Age_Weeks",
           "Last_cell_type", "Last_cell_type_num", "Last_cell_type_num2", "Cell_lineage", "TCRBCRlabel"]
for column in columns:
    result.obs[column] = annotations[column].to_numpy()
result.obs["Cellname"] = result.obs_names
result.obsm["X_umap"] = annotations[["UMAP1", "UMAP2"]].to_numpy(dtype=np.float32)
save_h5ad(result, GLOBAL / "Fetal_Immune_Atlas.h5ad")
