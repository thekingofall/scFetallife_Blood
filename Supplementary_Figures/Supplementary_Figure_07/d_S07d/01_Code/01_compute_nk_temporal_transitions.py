from pathlib import Path
import os
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT / "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_07/01_NK_Trajectories/Fetal_NK_CellRank_Input.h5ad"
adata = sc.read_h5ad(SOURCE)
import cellrank as cr
from cellrank.kernels import RealTimeKernel
from moscot.problems.time import TemporalProblem
import jax
import jax.numpy as jnp
if not hasattr(jnp, "DeviceArray"):
    jnp.DeviceArray = jax.Array
from scipy.sparse import csr_matrix, save_npz
scv.settings.set_figure_params("scvelo")
adata.uns["Last_cell_type_colors"] = ["#C71000FF", "#008EA0FF", "#8A4198FF"]
sc.settings.set_figure_params(frameon=False, dpi=100)
cr.settings.verbosity = 2
adata.obs["day"] = pd.to_numeric(adata.obs["day"].astype(str))
tp = TemporalProblem(adata)
tp = tp.prepare(time_key="day")
tp = tp.solve(epsilon=1e-3, tau_a=0.95, scale_cost="mean")
adata.obs["day"] = adata.obs["day"].astype(float).astype("category")
tmk = RealTimeKernel.from_moscot(tp)
tmk.compute_transition_matrix(
    self_transitions="all", conn_weight=0.2, threshold="auto"
)

flow_specs = (
    (
        "CX3CR1+ NK",
        adata.obs["Last_cell_type"].value_counts()[0:40].index.to_list(),
        "S07d_CX3CR1_terminal.pdf",
    ),
    (
        "CXCR6+ NK",
        adata.obs["Last_cell_type"].value_counts()[0:40].index.to_list(),
        "S07b_CXCR6_to_CX3CR1.pdf",
    ),
    (
        "CD56highCD16low NK",
        adata.obs["Last_cell_type"].unique(),
        "S07c_CD56high_to_CX3CR1.pdf",
    ),
)

for cluster, clusters, filename in flow_specs:
    ax = tmk.plot_single_flow(
        cluster_key="Last_cell_type",
        time_key="day",
        cluster=cluster,
        min_flow=0.5,
        xticks_step_size=4,
        show=False,
        clusters=clusters,
    )
    tick_ages = np.array([11.6, 21.7, 26.9, 32.4, 37.9])
    current_ages = np.array([float(label.get_text()) for label in ax.get_xticklabels()])
    tick_positions = np.interp(tick_ages, current_ages, ax.get_xticks())
    ax.set_xticks(tick_positions)
    ax.set_xticklabels([f"{age:.1f}" for age in tick_ages], rotation=90)
    ax.set_xlabel("pcw")
    ax.set_ylabel("Transition")
    letter = filename[3]
    output = ROOT / "Supplementary_Figures/Supplementary_Figure_07" / f"{letter}_S07{letter}" / "02_Figures"
    output.mkdir(parents=True, exist_ok=True)
    plt.savefig(output / filename, bbox_inches="tight")
    plt.savefig((output / filename).with_suffix(".png"), bbox_inches="tight", dpi=300)
    plt.close()

data_dir = ROOT / "Supplementary_Figures/Supplementary_Figure_07/b_S07b/03_PlotData"
save_npz(data_dir / "S07_NK_transition_matrix.npz", csr_matrix(tmk.transition_matrix))
adata.obs.to_csv(data_dir / "S07_NK_transition_cell_metadata.csv.gz")
