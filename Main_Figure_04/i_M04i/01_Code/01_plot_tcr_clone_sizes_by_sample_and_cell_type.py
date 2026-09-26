from pathlib import Path


SCF_PROJECT_ROOT = Path(__file__).resolve().parents[3]
SCF_GLOBAL = SCF_PROJECT_ROOT / "00_Global_Data"
SCF_SHARED = SCF_GLOBAL / "Shared_Inputs"
import matplotlib.pyplot as plt
import pandas as pd
from scirpy.pl import base


CLEAN_ROOT = SCF_PROJECT_ROOT
INPUT_DIR = SCF_SHARED / "Main_Figure_04_05/01_Chain_Pairing_Clonotypes"


def ordered_wide(data: pd.DataFrame, group_col: str) -> pd.DataFrame:
    data = data.copy()
    size_order = sorted(
        data["clone_size"].astype(str).unique(),
        key=lambda value: (value.startswith(">="), int(value.replace(">=", ""))),
    )
    data["clone_size"] = pd.Categorical(
        data["clone_size"].astype(str), categories=size_order, ordered=True
    )
    if "Main_Organ" in data.columns:
        organ_order = pd.CategoricalDtype(
            ["PBMC", "Liver", "Thymus", "Spleen"], ordered=True
        )
        data["Main_Organ"] = data["Main_Organ"].astype(organ_order)
        data = data.sort_values(["Main_Organ", group_col, "clone_size"])
    wide = data.pivot_table(
        index=group_col,
        columns="clone_size",
        values="proportion",
        aggfunc="sum",
        observed=False,
    ).fillna(0)
    return wide.loc[:, size_order]


def render(input_name: str, group_col: str, panel_id: str, folder: str) -> None:
    data = pd.read_csv(INPUT_DIR / input_name)
    wide = ordered_wide(data, group_col)

    output_dir = SCF_PROJECT_ROOT / "Main_Figure_04" / f"{panel_id[-1]}_{panel_id}" / "02_Figures"
    output_dir.mkdir(parents=True, exist_ok=True)


    ax = base.bar(wide, fig_kws={"figsize": (10, 10)})
    ax.figure.savefig(output_dir / f"{panel_id}.pdf", bbox_inches="tight")
    ax.figure.savefig(output_dir / f"{panel_id}.png", bbox_inches="tight", dpi=300)
    plt.close(ax.figure)
    data.to_csv(output_dir / f"{panel_id}_data.csv", index=False)


render(
    "M04i_TCR_clone_size_by_sample.csv",
    "MainID",
    "M04i_TCR_clone_size_by_sample",
    "M04i_TCR_clone_size_by_sample",
)
render(
    "M04j_TCR_clone_size_by_celltype.csv",
    "cell_type",
    "M04j_TCR_clone_size_by_celltype",
    "M04j_TCR_clone_size_by_celltype",
)
print("M04i and M04j complete")
