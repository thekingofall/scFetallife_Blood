from pathlib import Path

SCF_PROJECT_ROOT = Path(__file__).resolve().parents[3]
SCF_GLOBAL = SCF_PROJECT_ROOT / "00_Global_Data"
SCF_SHARED = SCF_GLOBAL / "Shared_Inputs"
import warnings

import matplotlib.pyplot as plt
import scanpy as sc


warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.filterwarnings("ignore", "No data for colormapping provided via", UserWarning)

clean_root = SCF_PROJECT_ROOT
m04h_dir = SCF_PROJECT_ROOT / "Main_Figure_04" / "h_M04h" / "02_Figures"
m05h_dir = SCF_PROJECT_ROOT / "Main_Figure_05" / "h_M05h" / "02_Figures"
m04h_dir.mkdir(parents=True, exist_ok=True)
m05h_dir.mkdir(parents=True, exist_ok=True)
m04h_data = m04h_dir.parent / "03_PlotData"
m05h_data = m05h_dir.parent / "03_PlotData"
m04h_data.mkdir(exist_ok=True)
m05h_data.mkdir(exist_ok=True)

plt.style.use("default")
plt.rcParams["font.family"] = "Arial"
plt.rcParams["figure.figsize"] = [4, 4]
plt.rcParams["figure.dpi"] = 300

adata = sc.read_h5ad(SCF_GLOBAL / "Fetal_Immune_Atlas.h5ad")


Tcell = adata[adata.obs.Cell_lineage.isin(["T/ILC"])]
Thy = Tcell[Tcell.obs.Main_Organ.isin(["Thymus"])]
DN = Thy[Thy.obs.Last_cell_type.isin(["DN(Q) T"])].copy()
DN.obs["Time_DN"] = (
    DN.obs["MainID"].astype("str")
    + "_"
    + DN.obs["Last_cell_type"].astype("str")
)
sc.pl.dotplot(DN, "DNTT", "Time_DN", swap_axes=True, show=False)
plt.savefig(m04h_dir / "M04h_DNTT_thymic_T.pdf", bbox_inches="tight")
plt.savefig(m04h_dir / "M04h_DNTT_thymic_T.png", bbox_inches="tight", dpi=300)
plt.close("all")
DN.obs[["MainID", "Last_cell_type", "Time_DN"]].to_csv(
    m04h_data / "M04h_DNTT_thymic_T_cell_census.csv"
)


Liver = adata[adata.obs.Main_Organ.isin(["Liver"])]
LiverB = Liver[Liver.obs.Last_cell_type.isin(["Pro-B"])].copy()
LiverB.obs["Time_ProB"] = (
    LiverB.obs["MainID"].astype("str")
    + "_"
    + LiverB.obs["Last_cell_type"].astype("str")
)
sc.pl.dotplot(LiverB, "DNTT", "Time_ProB", swap_axes=True, show=False)
plt.savefig(m05h_dir / "M05h_DNTT_liver_ProB.pdf", bbox_inches="tight")
plt.savefig(m05h_dir / "M05h_DNTT_liver_ProB.png", bbox_inches="tight", dpi=300)
plt.close("all")
LiverB.obs[["MainID", "Last_cell_type", "Time_ProB"]].to_csv(
    m05h_data / "M05h_DNTT_liver_ProB_cell_census.csv"
)

print("M04h and M05h complete")
