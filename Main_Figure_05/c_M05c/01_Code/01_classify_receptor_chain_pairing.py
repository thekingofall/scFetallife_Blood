from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scirpy as ir

root = Path(__file__).resolve().parents[3]
panel = Path(__file__).resolve().parents[1]
data_dir = panel / "03_PlotData"
data_dir.mkdir(parents=True, exist_ok=True)
receptor = "BCR"
atlas = ad.read_h5ad(root / "00_Global_Data/Fetal_Immune_Atlas.h5ad", backed="r")
obs = atlas.obs.copy()
umap = np.asarray(atlas.obsm["X_umap"])
atlas.file.close()
contigs = pd.read_csv(root / f"00_Global_Data/Shared_Inputs/Main_Figure_04_05/00_Receptor_Annotations/01_Aggregated_Contigs/Gao_All_{receptor}.csv", index_col=0)
contigs["barcode"] = contigs["barcodename"]
contigs["productive"] = contigs["productive"].astype(str)
aligned = data_dir / f"{receptor}_contig_annotations.csv"
contigs.to_csv(aligned, index=False)
vdj = ir.io.read_10x_vdj(aligned)
master = ad.AnnData(X=sp.csr_matrix((len(obs), 0)), obs=obs)
master.obsm["X_umap"] = umap
ir.pp.merge_with_ir(master, vdj)
ir.tl.chain_qc(master)
result = master.obs[["MainID", "Main_Organ", "Last_cell_type", "chain_pairing", "receptor_subtype"]].copy()
result["UMAP1"] = umap[:, 0]
result["UMAP2"] = umap[:, 1]
result = result[result["chain_pairing"].astype(str) != "no IR"]
result.index.name = "Cellname"
result.to_csv(data_dir / f"{receptor}_chain_pairing_UMAP.csv.gz")
