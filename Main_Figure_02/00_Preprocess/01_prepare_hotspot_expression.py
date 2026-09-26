"""Align raw UMI counts with the annotated atlas and prepare early/late PBMC inputs."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "00_Settings"))

import numpy as np
import pandas as pd
import scanpy as sc
from preprocess_inputs import GLOBAL, read_counts, save_h5ad

here = Path(__file__).resolve().parent
atlas = sc.read_h5ad(GLOBAL / "Fetal_Immune_Atlas.h5ad")
genes = pd.Index(pd.read_csv(here / "01_Inputs/Hotspot_Genes.csv")["gene"])
genes = genes[~genes.str.startswith(("MT", "RP"))]
data = atlas[atlas.obs["Main_Organ"].astype(str).eq("PBMC"), genes].copy()
counts = read_counts(data.obs_names, data.var_names)
data.layers["counts"] = counts.X.copy()
data.layers["counts_csc"] = counts.X.tocsc()
data.layers["log_normalized"] = data.X.copy()
data.obs["Stage"] = np.where(data.obs["Post_Conception_Age_Weeks"] <= 26, "Early", "Late")
for stage in ["Early", "Late"]:
    subset = data[data.obs["Stage"] == stage].copy()
    save_h5ad(subset, here / "02_Processed" / f"M02_PBMC_{stage}_Hotspot.h5ad")
