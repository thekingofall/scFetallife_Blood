from pathlib import Path


SCF_PROJECT_ROOT = Path(__file__).resolve().parents[3]
SCF_GLOBAL = SCF_PROJECT_ROOT / "00_Global_Data"
SCF_SHARED = SCF_GLOBAL / "Shared_Inputs"
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


CLEAN_ROOT = SCF_PROJECT_ROOT
INPUT_DIR = SCF_SHARED / "Main_Figure_04_05/01_Chain_Pairing_Clonotypes"
PALETTE = ["#C71000FF", "#f49128", "#023f75", "#5A9599FF"]
ORGAN_ORDER = ["PBMC", "Liver", "Thymus", "Spleen"]


def render(input_name: str, panel_id: str, figure_number: str, figsize: tuple[float, float], dpi: int) -> None:
    source = pd.read_csv(INPUT_DIR / input_name)
    source["Main_Organ"] = pd.Categorical(
        source["Main_Organ"], categories=ORGAN_ORDER, ordered=True
    )
    plot_adata = ad.AnnData(
        X=np.zeros((source.shape[0], 1), dtype=np.float32),
        obs=source.set_index("cell_id")[["Main_Organ"]].copy(),
    )
    plot_adata.obsm["X_umap"] = source[["UMAP1", "UMAP2"]].to_numpy()

    output_dir = (
        SCF_PROJECT_ROOT
        / f"Main_Figure_{figure_number}"
        / f"{panel_id[-1]}_{panel_id}"
        / "02_Figures"
    )
    output_dir.mkdir(parents=True, exist_ok=True)


    plt.style.use("default")
    plt.rcParams["figure.figsize"] = list(figsize)
    plt.rcParams["figure.dpi"] = dpi
    sc.pl.umap(
        plot_adata,
        color="Main_Organ",
        palette=PALETTE,
        show=False,
    )
    plt.savefig(output_dir / f"{panel_id}_organ_UMAP.pdf", bbox_inches="tight")
    plt.savefig(
        output_dir / f"{panel_id}_organ_UMAP.png",
        bbox_inches="tight",
        dpi=300,
    )
    plt.close("all")


render(
    "M04c_TCR_Chain_Pairing_UMAP.csv",
    "M04d",
    "04",
    (2, 2),
    250,
)
render(
    "M05c_BCR_Chain_Pairing_UMAP.csv",
    "M05d",
    "05",
    (4, 4),
    300,
)
print("M04d and M05d complete")
