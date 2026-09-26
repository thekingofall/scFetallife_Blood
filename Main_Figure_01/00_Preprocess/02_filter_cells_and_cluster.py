"""Remove doublets, apply RNA quality filters and calculate PCA, UMAP and Leiden clusters."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "00_Settings"))

import numpy as np
import pandas as pd
import scanpy as sc
from preprocess_inputs import GLOBAL, sample_inputs, save_h5ad

output = Path(__file__).resolve().parent / "02_Processed"
data = sc.read_h5ad(output / "M01_CellTypist.h5ad")
metadata = sample_inputs().set_index("Name")
for column in ["MainID", "Main_Organ", "Gestational_Age_Weeks", "Post_Conception_Age_Weeks"]:
    data.obs[column] = data.obs["Name"].map(metadata[column])
contigs = GLOBAL / "Shared_Inputs/Main_Figure_04_05/00_Receptor_Annotations/01_Aggregated_Contigs"
tcr = pd.read_csv(contigs / "Gao_All_TCR.csv", usecols=["barcodename"])
bcr = pd.read_csv(contigs / "Gao_All_BCR.csv", usecols=["barcodename"])
t = data.obs_names.isin(tcr["barcodename"])
b = data.obs_names.isin(bcr["barcodename"])
data.obs["TCRBCRlabel"] = np.select([t & b, t, b], ["Both", "TCR", "BCR"], default="None")
data = data[(data.obs["TCRBCRlabel"] != "Both") & ~data.obs["predicted_doublets"].astype(bool)].copy()
sc.pp.filter_cells(data, min_genes=200)
sc.pp.filter_genes(data, min_cells=3)
data.var["mt"] = data.var_names.str.startswith("MT-")
# Calculate QC metrics on the library-normalized, log-transformed matrix.
sc.pp.calculate_qc_metrics(data, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
data = data[(data.obs["n_genes_by_counts"] < 8000) & (data.obs["pct_counts_mt"] < 10)].copy()
data.uns["log1p"]["base"] = None
data.raw = data.copy()
sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5)
data = data[:, data.var["highly_variable"]].copy()
sc.pp.regress_out(data, ["total_counts", "pct_counts_mt"])
sc.pp.scale(data, max_value=10)
sc.tl.pca(data, svd_solver="arpack")
sc.pp.neighbors(data, n_neighbors=10, n_pcs=40)
sc.tl.umap(data)
sc.tl.leiden(data)
save_h5ad(data, output / "M01_Filtered_Clustered.h5ad")
