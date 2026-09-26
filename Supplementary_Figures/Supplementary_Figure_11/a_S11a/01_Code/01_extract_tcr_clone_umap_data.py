"""Export TCR clone annotations and the stored UMAP coordinates."""
from pathlib import Path
import anndata as ad
import pandas as pd

root=Path(__file__).resolve().parents[4]
adata=ad.read_h5ad(root/"00_Global_Data/Fetal_Immune_Atlas_TCR.h5ad",backed="r")
try:
    fields=["clone_id","clone_id_size","clonal_expansion","MainID","Main_Organ","Last_cell_type"]
    output=adata.obs[fields].copy()
    output.insert(0,"cell_id",output.index.astype(str))
    output["UMAP1"]=adata.obsm["X_umap"][:,0]
    output["UMAP2"]=adata.obsm["X_umap"][:,1]
finally:
    adata.file.close()
dest=Path(__file__).resolve().parents[1]/"03_PlotData"
dest.mkdir(parents=True,exist_ok=True)
output.to_csv(dest/"S11a_TCR_clone_annotations.csv.gz",index=False)
