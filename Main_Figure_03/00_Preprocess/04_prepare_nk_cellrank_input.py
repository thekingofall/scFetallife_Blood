"""Prepare the fitted NK velocity object for the temporal CellRank analysis."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "00_Settings"))

import anndata as ad
import h5py
from anndata._io.specs import read_elem
import pandas as pd
from preprocess_inputs import GLOBAL, save_h5ad
atlas = ad.read_h5ad(GLOBAL / "Fetal_Immune_Atlas.h5ad", backed="r")
with h5py.File(GLOBAL / "Raw_Inputs/NK/Cellrank_tpNK.h5ad") as h:
    data = ad.AnnData(X=read_elem(h["X"]), obs=read_elem(h["obs"]), var=read_elem(h["var"]),
        obsm=read_elem(h["obsm"]), obsp=read_elem(h["obsp"]), uns=read_elem(h["uns"]),
        layers={"velocity": read_elem(h["layers/velocity"])})
data = data[data.obs_names.isin(atlas.obs_names)].copy()
data.obs = atlas.obs.loc[data.obs_names].copy()
data.obs["day"] = pd.Categorical(data.obs["Post_Conception_Age_Weeks"])
save_h5ad(data, GLOBAL / "Shared_Inputs/Supplementary_Figures/Supplementary_Figure_07/01_NK_Trajectories/Fetal_NK_CellRank_Input.h5ad")
atlas.file.close()
