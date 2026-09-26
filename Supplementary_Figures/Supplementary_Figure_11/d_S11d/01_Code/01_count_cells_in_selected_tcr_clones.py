"""Tabulate cells belonging to the 15 manuscript TCR clonotypes."""
from pathlib import Path
import anndata as ad
import pandas as pd

root=Path(__file__).resolve().parents[4]
panel=Path(__file__).resolve().parents[1]
selected=pd.read_csv(root/"00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_11/d_S11d/S11d_Top15_Clones_by_Sample_and_Cell_Type.csv",dtype={"clone_id":str})
obs=pd.read_csv(root/"Supplementary_Figures/Supplementary_Figure_11/b_S11b/03_PlotData/TCR_clonotype_network_coordinates.csv.gz",dtype={"clone_id":str})
adata=ad.read_h5ad(root/"00_Global_Data/Fetal_Immune_Atlas_TCR.h5ad",backed="r")
try:
    obs["Last_cell_type"]=obs.cell_id.map(adata.obs.Last_cell_type.astype(str))
finally:
    adata.file.close()
obs["clone_id"]=obs.clone_id.astype(str)
obs=obs[obs.clone_id.isin(selected.clone_id.unique())]
table=obs.groupby(["clone_id","MainID","Last_cell_type"],observed=True).size().reset_index(name="cells")
order={v:i for i,v in enumerate(selected.clone_id.unique())}
table=table.sort_values("clone_id",key=lambda s:s.map(order))
out=panel/"03_PlotData";out.mkdir(parents=True,exist_ok=True)
table.to_csv(out/"S11d_TCR_top15_plot_source_data.csv",index=False)
