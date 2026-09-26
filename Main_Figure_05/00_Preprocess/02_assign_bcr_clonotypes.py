"""Infer BCR clones using the heavy-chain distance threshold."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "00_Settings"))
import numpy as np
import scanpy as sc
import dandelion as ddl
from preprocess_inputs import GLOBAL, save_h5ad
output = Path(__file__).resolve().parent / "02_Processed"
output.mkdir(exist_ok=True)
data = sc.read_h5ad(GLOBAL / "Fetal_Immune_Atlas_BCR.h5ad")
receptors = ddl.from_scirpy(data)
for column in ["v_call", "j_call", "junction_aa"]:
    receptors.data[column] = receptors.data[column].replace("", np.nan)
receptors.data = receptors.data.dropna(subset=["v_call", "j_call", "junction_aa"])
receptors.data["junction_length"] = receptors.data["junction_aa"].str.len()
receptors.threshold = 0.15858785672187678
ddl.tl.define_clones(receptors, key_added="changeo_clone_id", model="hh_s5f")
ddl.tl.generate_network(receptors, key="sequence_alignment", layout_method="sfdp", clone_key="changeo_clone_id")
ddl.tl.clone_size(receptors, clone_key="changeo_clone_id", max_size=6)
ddl.tl.transfer(data, receptors, clone_key="changeo_clone_id")
receptors.data.to_csv(output / "M05_BCR_Clone_Annotations.csv", index=False)
receptors.write_pkl(str(output / "M05_BCR_Clonotypes.pkl"))
save_h5ad(data, output / "M05_BCR_Clonotypes.h5ad")
