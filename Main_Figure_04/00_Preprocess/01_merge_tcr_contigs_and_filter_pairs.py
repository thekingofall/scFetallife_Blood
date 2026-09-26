"""Attach TCR contigs to atlas cells and retain paired receptors for clone analysis."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "00_Settings"))

import pandas as pd
import scanpy as sc
import scirpy as ir
from preprocess_inputs import GLOBAL, save_h5ad
here = Path(__file__).resolve().parent
output = here / "02_Processed"
output.mkdir(exist_ok=True)
data = sc.read_h5ad(GLOBAL / "Fetal_Immune_Atlas.h5ad")
contigs = pd.read_csv(GLOBAL / "Shared_Inputs/Main_Figure_04_05/00_Receptor_Annotations/01_Aggregated_Contigs/Gao_All_TCR.csv", low_memory=False)
contigs = contigs[contigs["barcodename"].isin(data.obs_names)].copy()
contigs["barcode"] = contigs["barcodename"]
contigs.to_csv(output / "M04_TCR_Contigs.csv", index=False)
receptors = ir.io.read_10x_vdj(output / "M04_TCR_Contigs.csv")
ir.pp.merge_with_ir(data, receptors)
ir.tl.chain_qc(data)
data = data[~data.obs["chain_pairing"].isin(["no IR", "multichain", "orphan VDJ", "orphan VJ"])].copy()
cells = pd.Index(pd.read_csv(here / "01_Inputs/TCR_Cells.csv.gz")["Cellname"])
if len(cells.difference(data.obs_names)):
    raise ValueError("Curated TCR cells are missing after chain pairing; check the contig input")
save_h5ad(data[cells].copy(), output / "M04_Paired_TCR.h5ad")
