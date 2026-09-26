from pathlib import Path
import os
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT / "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_07/01_NK_Trajectories/Fetal_NK_CellRank_Input.h5ad"
adata = sc.read_h5ad(SOURCE)
OUTPUT = ROOT / "Supplementary_Figures/Supplementary_Figure_07/a_S07a/02_Figures"
OUTPUT.mkdir(parents=True, exist_ok=True)
scv.settings.set_figure_params("scvelo")
# scVelo 0.2.x uses a pandas categorical setter removed in pandas 2.
import importlib
import inspect
legend_utils = importlib.import_module("scvelo.plotting.utils")
legend_source = inspect.getsource(legend_utils.set_legend)
legacy_setter = "obs_vals.cat.categories = obs_vals.cat.categories.astype(str)"
if legacy_setter in legend_source:
    exec(legend_source.replace(legacy_setter, "obs_vals = obs_vals.cat.rename_categories(obs_vals.cat.categories.astype(str))"), legend_utils.__dict__)
    importlib.import_module("scvelo.plotting.scatter").set_legend = legend_utils.set_legend

col1 = ["#C71000FF", "#008EA0FF", "#8A4198FF"]
scv.pl.velocity_embedding_stream(
    adata,
    basis="X_umap",
    color="Last_cell_type",
    palette=col1,
    size=20,
    alpha=0.8,
    legend_loc="right margin",
    figsize=(7, 5),
    legend_fontsize=9,
    show=False,
    title="",
)
png = OUTPUT / "S07a_NK_RNA_velocity.png"
plt.savefig(png, bbox_inches="tight", dpi=300)
try:
    plt.savefig(OUTPUT / "S07a_NK_RNA_velocity.pdf", bbox_inches="tight")
except ValueError:
    from PIL import Image
    Image.open(png).convert("RGB").save(OUTPUT / "S07a_NK_RNA_velocity.pdf", resolution=300.0)
plt.close()
table = adata.obs.copy()
table[["UMAP1", "UMAP2"]] = adata.obsm["X_umap"]
table[["velocity_UMAP1", "velocity_UMAP2"]] = adata.obsm["velocity_umap"]
table.to_csv(OUTPUT.parent / "03_PlotData/S07a_NK_velocity_coordinates.csv.gz")
