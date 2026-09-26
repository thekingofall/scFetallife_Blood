"""Align the spliced/unspliced NK-cell object with the atlas cell annotations."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "00_Settings"))

import anndata as ad
from preprocess_inputs import GLOBAL, save_h5ad
atlas = ad.read_h5ad(GLOBAL / "Fetal_Immune_Atlas.h5ad", backed="r")
data = ad.read_h5ad(GLOBAL / "Raw_Inputs/NK/SubNK.h5ad")
if not {"spliced", "unspliced"}.issubset(data.layers):
    raise ValueError("The NK source must contain spliced and unspliced count layers")
data = data[data.obs_names.isin(atlas.obs_names)].copy()
data.obs = atlas.obs.loc[data.obs_names].copy()
save_h5ad(data, GLOBAL / "Shared_Inputs/Supplementary_Figures/Supplementary_Figure_07/01_NK_Trajectories/Fetal_NK_Velocity.h5ad")
atlas.file.close()
