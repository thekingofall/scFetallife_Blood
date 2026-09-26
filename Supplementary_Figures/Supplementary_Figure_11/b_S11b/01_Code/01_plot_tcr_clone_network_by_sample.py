"""Plot exact TCR clonotype networks with Scirpy 0.12."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir

ORGAN_ORDER = ["PBMC", "Liver", "Thymus", "Spleen"]

ORGAN_COLORS = {"PBMC": "#C71000", "Liver": "#f49128", "Thymus": "#023f75", "Spleen": "#5A9599"}

SAMPLE_COLORS = {
    "B11.6_P24": "#ccff99", "B17.4_P3": "#868686", "B18.0_P4": "#f6313e", "B18.6_P5": "#ffc179", "B20.9_P7": "#001588", "B21.7_P8": "#cccc33", "B22.4_P9": "#993333", "B22.4_P23": "#006d2c", "B23.4_P10": "#b5ad64", "B24.6_P11": "#8f1336", "B26.9_P12": "#ffa300", "B29.1_P13": "#333399", "B31.3_P15": "#98d9e9", "B32.4_P16": "#cc9966", "B33.3_P17": "#1f78b4", "B34.1_P18": "#ba7fd0", "B36.1_P19": "#009966", "B37.9_P20": "#ff6666", "B38.1_P21": "#b2df8a", "B39.1_P22": "#9da8e2", "L10.1_P25": "#999933", "L18.6_P5": "#a6cee3", "L24.6_P11": "#cc0033", "T10.0_P1": "#663366", "T10.1_P25": "#00af99", "T18.6_P5": "#46a040", "T24.6_P11": "#ff5a00", "S18.6_P5": "#333366", "S24.6_P11": "#de77ae",
}

PANEL_SPECS = {
    "clone_id": {"stem": "S29A_clone_ID", "inches": (14.5, 8.0)},
    "clone_organ": {"stem": "S29A_clone_Organ", "inches": (9.5, 7.0)},
    "organ_count_a": {"stem": "S29A_organ_count", "inches": (10.6667, 4.0)},
    "organ_count_c": {"stem": "S29C_organ_count", "inches": (10.6667, 4.0)},
}

def restore_numeric_clone_index_order(data) -> None:
    """Order clonotype indices numerically before calculating the Scirpy network."""
    cell_indices = data.uns["clone_id"]["cell_indices"]
    count = len(cell_indices)
    expected = {str(index) for index in range(count)}
    if set(cell_indices) != expected:
        raise ValueError("clone_id cell_indices keys are not a complete numeric range")
    data.uns["clone_id"]["cell_indices"] = {str(index): cell_indices[str(index)] for index in range(count)}

def _plot_network(data, color: str, panel_key: str, pdf: Path, png: Path) -> None:
    import scirpy as ir

    width, height = PANEL_SPECS[panel_key]["inches"]
    if color == "Main_Organ":
        palette = [ORGAN_COLORS[value] for value in ORGAN_ORDER]
        data.obs[color] = pd.Categorical(data.obs[color].astype(str), categories=ORGAN_ORDER, ordered=True)
        panel_size = (7.0, 7.0)
        label_fontsize = 5
    else:
        sample_order = sorted(data.obs[color].astype(str).unique())
        if set(sample_order) != set(SAMPLE_COLORS):
            raise ValueError("Sample IDs differ from the sample palette")
        data.obs[color] = pd.Categorical(data.obs[color].astype(str), categories=sample_order, ordered=True)
        palette = [SAMPLE_COLORS[value] for value in sample_order]
        panel_size = (12.0, 8.0)
        label_fontsize = 6
    ax = ir.pl.clonotype_network(
        data,
        color=color,
        base_size=20,
        label_fontsize=label_fontsize,
        panel_size=panel_size,
        palette=palette,
    )
    figure = ax.figure
    figure.set_size_inches(width, height, forward=True)
    figure.savefig(pdf, format="pdf", facecolor="white")
    figure.savefig(png, format="png", dpi=300, facecolor="white")
    plt.close(figure)

root=Path(__file__).resolve().parents[4]
panel=Path(__file__).resolve().parents[1]
data=sc.read_h5ad(root/"00_Global_Data/Fetal_Immune_Atlas_TCR.h5ad")
ir.pp.ir_dist(data,sequence="aa")
ir.tl.define_clonotypes(data,receptor_arms="all",dual_ir="primary_only")
ir.tl.clonal_expansion(data)
restore_numeric_clone_index_order(data)
ir.tl.clonotype_network(data,min_cells=2,random_state=42)
coords=np.asarray(data.obsm["X_clonotype_network"])
selected=~pd.isna(coords).all(axis=1)
table=data.obs.loc[selected,["clone_id","clone_id_size","MainID","Main_Organ"]].copy()
table.insert(0,"cell_id",table.index.astype(str))
table["x"]=coords[selected,0];table["y"]=coords[selected,1]
output=panel/"03_PlotData";output.mkdir(parents=True,exist_ok=True)
table.to_csv(output/"TCR_clonotype_network_coordinates.csv.gz",index=False)
figures=panel/"02_Figures";figures.mkdir(parents=True,exist_ok=True)
key="clone_id" if panel.name=="b_S11b" else "clone_organ"
color="MainID" if key=="clone_id" else "Main_Organ"
stem="S11b_TCR_clone_network_by_sample" if key=="clone_id" else "S11c_TCR_clone_network_by_organ"
plt.rcParams["font.family"]="Arial"
_plot_network(data,color,key,figures/(stem+".pdf"),figures/(stem+".png"))
print(f"{len(table)} cells in {table.clone_id.nunique()} clonotypes")
