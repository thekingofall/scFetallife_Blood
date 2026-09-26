from pathlib import Path


SCF_PROJECT_ROOT = Path(__file__).resolve().parents[4]
SCF_GLOBAL = SCF_PROJECT_ROOT / "00_Global_Data"
SCF_SHARED = SCF_GLOBAL / "Shared_Inputs"
import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd
import scanpy as sc


CLEAN_ROOT = SCF_PROJECT_ROOT
OUTPUT_DIR = SCF_PROJECT_ROOT / "Supplementary_Figures/Supplementary_Figure_01" / "a_S01a" / "02_Figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams["font.family"] = "Arial"

adata1 = sc.read_h5ad(SCF_GLOBAL / "Fetal_Immune_Atlas.h5ad")

categories = pd.Categorical(adata1.obs["Last_cell_type_num"])
sorted_categories = sorted(categories.categories, key=lambda x: int(x.split("_")[0]))
adata1.obs["Last_cell_type_num"] = pd.Categorical(
    adata1.obs["Last_cell_type_num"],
    categories=sorted_categories,
    ordered=True,
)
adata1.obs["Cell_type"] = adata1.obs["Last_cell_type_num"].cat.rename_categories(
    lambda label: label.split("_", 1)[1]
)


Mainmarker = [
    "PTPRC", "PLAC8", "SERPINB1", "HLA-DRA", "KIT", "GATA2", "CD34", "SPINK2",
    "JCHAIN", "IGLL1", "CD79B", "TCL1A", "IGKC", "MS4A1", "CD19", "CXCR5",
    "DNTT", "CDC45", "DHFR", "CD27", "CD5", "TNFRSF17", "LTB", "CD3E", "CD7",
    "IL32", "CCR7", "CD8A", "CD4", "NKG7", "KLRB1", "KLRD1", "KLRC1", "XCL2",
    "CXCR6", "CX3CR1", "GZMH", "GNLY", "GZMB", "GZMK", "NCAM1", "TRDC", "TRDV1",
    "TRDV2", "TRGV9", "FOXP3", "PDCD1", "GNG4", "CCR9", "MPO", "CCL4", "FCGR3A",
    "LYZ", "S100A9", "CD14", "IL3RA", "CLEC9A", "CD1C", "C1QA", "VCAM1", "TPSAB1",
    "PF4", "ITGA2B", "ESAM", "CD177", "GYPA", "GATA1", "KLF1", "ALAS2", "HBA1",
    "BPGM", "HBE1", "MKI67",
]


colname4 = [
    "#C71000FF", "#008EA0FF", "#8A4198FF", "#5A9599FF", "#FF6348FF", "#84D7E1FF",
    "#FF95A8FF", "#3D3B25FF", "#ADE2D0FF", "#1A5354FF", "#3F4041FF", "#fa6e01",
    "#972b1d", "#e6a84b", "#4c211b", "#ff717f", "#009966", "#c62d17", "#023f75",
    "#ea894e", "#266b69", "#eb4601", "#f6c619", "#f49128", "#194a55", "#c29f62",
    "#83ba9e", "#187c65", "#A6CEE3", "#223e9c", "#aebea6", "#edae11", "#c74732",
    "#6a73cf", "#edd064", "#0eb0c8", "#f2ccac", "#a1d5b9", "#e1abbc", "#46A040",
    "#00AF99", "#FFC179", "#98D9E9", "#F6313E", "#FFA300", "#333366", "#663366",
    "#FF6666", "#8F1336", "#0081C9", "#CC0033", "#CC9966", "#CC0033", "#999933",
    "#CCCC33", "#CCFF99", "#333399", "#001588", "#490C65", "#BA7FD0", "#1F78B4",
    "#DE77AE", "#2f2f2f", "#006D2C", "#868686", "#9DA8E2", "#91C392", "#FF9900",
    "#339966", "#993333",
]


custom_fontsize = 25
custom_rc = {
    "font.family": "Arial",
    "font.size": custom_fontsize,
    "axes.titlesize": custom_fontsize,
    "axes.labelsize": custom_fontsize,
    "xtick.labelsize": custom_fontsize,
    "ytick.labelsize": custom_fontsize,
}
with plt.rc_context(rc=custom_rc):
    sc.pl.stacked_violin(
        adata1,
        Mainmarker,
        "Cell_type",
        show=False,
        swap_axes=False,


        row_palette=colname4,
        standard_scale="var",
    )
plt.savefig(OUTPUT_DIR / "S01a_marker_stacked_violin.pdf")


plt.savefig(OUTPUT_DIR / "S01a_marker_stacked_violin.png", dpi=150)
plt.close("all")

pd.DataFrame({"marker": Mainmarker}).to_csv(
    SCF_PROJECT_ROOT / "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_01/a_S01a/03_PlotData/S01a_marker_order.csv", index=False
)
adata1.obs[["MainID", "Last_cell_type_num"]].value_counts().rename("cells").reset_index().to_csv(
    OUTPUT_DIR / "S01a_cell_census.csv", index=False
)
print("S01a complete")
