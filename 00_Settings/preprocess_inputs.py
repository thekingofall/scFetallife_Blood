"""Shared inputs for the figure preprocessing scripts."""
from pathlib import Path
import anndata as ad
import pandas as pd
import scanpy as sc
ROOT = Path(__file__).resolve().parents[1]
GLOBAL = ROOT / "00_Global_Data"
INPUTS = ROOT / "Main_Figure_01/00_Preprocess/01_Inputs"

def sample_inputs():
    return pd.read_csv(INPUTS / "Sample_Inputs.csv")

def read_sample_counts(row):
    path = ROOT / row.matrix_file
    if not path.is_file():
        raise FileNotFoundError(f"Place the Cell Ranger matrix at {path}")
    data = sc.read_10x_h5(path)
    data.var_names_make_unique()
    data.obs_names = str(row.Name) + "_" + data.obs_names
    data.obs["Name"] = str(row.Name)
    return data

def read_counts(cells, genes):
    blocks = []
    for row in sample_inputs().itertuples(index=False):
        data = read_sample_counts(row)
        data = data[data.obs_names.isin(cells)].copy()
        if data.n_obs:
            blocks.append(data)
    merged = ad.concat(blocks, join="outer", merge="first", fill_value=0)
    missing_cells = cells.difference(merged.obs_names)
    missing_genes = genes.difference(merged.var_names)
    if len(missing_cells) or len(missing_genes):
        raise ValueError(f"Raw matrices lack {len(missing_cells)} cells and {len(missing_genes)} genes")
    return merged[cells, genes].copy()

def save_h5ad(data, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.stem + ".partial.h5ad")
    data.write_h5ad(temporary, compression="gzip")
    temporary.replace(path)
