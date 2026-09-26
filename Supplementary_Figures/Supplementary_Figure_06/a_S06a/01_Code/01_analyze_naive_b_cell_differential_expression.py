"""Differential expression between naive B-cell subsets and across fetal tissues."""
from __future__ import annotations

import argparse
from pathlib import Path


SCF_PROJECT_ROOT = Path(__file__).resolve().parents[4]
SCF_GLOBAL = SCF_PROJECT_ROOT / "00_Global_Data"
SCF_SHARED = SCF_GLOBAL / "Shared_Inputs"
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


CANONICAL_ATLAS = SCF_GLOBAL / "Fetal_Immune_Atlas.h5ad"
DEFAULT_OUTPUT = SCF_PROJECT_ROOT / "Supplementary_Figures/Supplementary_Figure_06"
CELL_TYPES = ("CXCR5+ Naïve B", "CXCR5- Naïve B")
DISPLAY_ORGANS = ("Liver", "Spleen", "PBMC")
RANKING_ORGANS = ("Liver", "PBMC", "Spleen", "Thymus")


def dense(matrix) -> np.ndarray:
    return np.asarray(matrix.toarray() if sparse.issparse(matrix) else matrix, dtype=float)


def validate_inputs(adata: ad.AnnData) -> pd.DataFrame:
    required = {"Last_cell_type", "Main_Organ", "MainID", "Name"}
    missing = required - set(adata.obs.columns)
    if missing:
        raise ValueError(f"Atlas is missing obs columns: {sorted(missing)}")
    counts = adata.obs["Last_cell_type"].astype(str).value_counts()
    if any(counts.get(cell_type, 0) == 0 for cell_type in CELL_TYPES):
        raise ValueError("Both naive B-cell subsets are required")
    return counts.reindex(CELL_TYPES).rename_axis("cell_type").reset_index(name="cell_count")



def rank_all_groups(subset: ad.AnnData, groupby: str) -> pd.DataFrame:
    subset = subset.copy()
    subset.uns.setdefault("log1p", {})["base"] = None
    sc.tl.rank_genes_groups(subset, groupby=groupby, method="wilcoxon", use_raw=False, n_genes=subset.n_vars)
    result = sc.get.rank_genes_groups_df(subset, group=None)
    result["group"] = result["group"].astype(str)
    return result


def plot_heatmap(matrix: pd.DataFrame, output_pdf: Path, output_png: Path, *, vmin: float, vmax: float,
                 colorbar_label: str, title: str, annotations: list[tuple[int, int, str]]) -> None:
    n_rows, n_cols = matrix.shape
    height = max(12.0, min(27.0, n_rows * 0.17))
    fig, ax = plt.subplots(figsize=(10.0, height))
    image = ax.imshow(matrix.to_numpy(), cmap="RdBu_r", vmin=vmin, vmax=vmax, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(n_cols), matrix.columns, rotation=90, fontsize=10)
    ax.set_yticks(np.arange(n_rows), matrix.index, fontsize=max(4.5, min(8.0, 470 / max(n_rows, 1))))
    ax.set_title(title, fontsize=15, pad=14)
    ax.tick_params(length=2.5, width=0.6)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    for start, end, label in annotations:
        midpoint = (start + end) / 2
        ax.plot([n_cols - 0.42, n_cols - 0.26, n_cols - 0.26, n_cols - 0.42],
                [start - 0.45, start - 0.45, end + 0.45, end + 0.45], color="black", lw=0.8, clip_on=False)
        ax.text(n_cols - 0.08, midpoint, label, va="center", ha="left", fontsize=8.5, clip_on=False)
    fig.subplots_adjust(left=0.20, right=0.50, bottom=0.08, top=0.95)
    colorbar_ax = fig.add_axes([0.88, 0.43, 0.018, 0.14])
    colorbar = fig.colorbar(image, cax=colorbar_ax)
    colorbar.set_label(colorbar_label, fontsize=9)
    fig.savefig(output_pdf)
    fig.savefig(output_png, dpi=300)
    plt.close(fig)


def build_s6a(master: ad.AnnData, panel: Path) -> None:
    subset = master[master.obs["Last_cell_type"].astype(str).isin(CELL_TYPES)].copy()
    subset.obs["naive_b_group"] = pd.Categorical(
        subset.obs["Last_cell_type"].astype(str), categories=CELL_TYPES, ordered=True
    )
    result = rank_all_groups(subset, "naive_b_group")
    result.to_csv(panel / "03_Source_Data" / "S6a_CXCR5_NaiveB_rank_genes_groups.csv", index=False)

    filtered = result.loc[
        result["logfoldchanges"].gt(3.0) & result["pvals_adj"].lt(0.05)
    ].head(100).copy()
    filtered.to_csv(
        SCF_PROJECT_ROOT / "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_06/a_S06a/03_PlotData/S6a_CXCR5_NaiveB_original_filter_top100.csv",
        index=False,
    )
    selected = filtered["names"].dropna().astype(str).tolist()
    selected = list(dict.fromkeys(gene for gene in selected if gene in subset.var_names))
    lookup = result.pivot_table(index="names", columns="group", values="logfoldchanges", aggfunc="first")
    matrix = lookup.reindex(selected).reindex(columns=CELL_TYPES).fillna(0.0)
    matrix.to_csv(SCF_PROJECT_ROOT / "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_06/a_S06a/03_PlotData/S6a_CXCR5_NaiveB_logfoldchange_matrix.csv")
    first_group = str(filtered.iloc[0]["group"])
    split = int(filtered["group"].astype(str).eq(first_group).cumprod().sum())
    plot_heatmap(
        matrix,
        panel / "02_Figures" / "FigureS6a_CXCR5_NaiveB_Differential_Expression.pdf",
        panel / "02_Figures" / "FigureS6a_CXCR5_NaiveB_Differential_Expression.png",
        vmin=-4.0,
        vmax=4.0,
        colorbar_label="Log fold change",
        title="CXCR5+ and CXCR5− naive B cells",
        annotations=[
            (0, split - 1, "Regulation of lymphocyte activation\nCytokine signalling in immune system"),
            (split, len(matrix) - 1, "MAPK family signalling cascades\nRegulation of developmental growth"),
        ],
    )


def build_s6b(master: ad.AnnData, panel: Path) -> None:
    ranking_subset = master[
        master.obs["Last_cell_type"].astype(str).isin(CELL_TYPES)
        & master.obs["Main_Organ"].astype(str).isin(RANKING_ORGANS)
    ].copy()
    ranking_subset.obs["organ_group"] = pd.Categorical(
        ranking_subset.obs["Main_Organ"].astype(str), categories=RANKING_ORGANS, ordered=True
    )
    if set(ranking_subset.obs["organ_group"].dropna().astype(str)) != set(RANKING_ORGANS):
        raise ValueError("S6b four-organ ranking background is incomplete")
    result = rank_all_groups(ranking_subset, "organ_group")
    result.to_csv(panel / "03_Source_Data" / "S6b_All_NaiveB_organ_rank_genes_groups.csv", index=False)
    genes = []
    for organ in DISPLAY_ORGANS:
        genes.extend(result.loc[result["group"].eq(organ), "names"].dropna().astype(str).head(20))
    genes = list(dict.fromkeys(gene for gene in genes if gene in ranking_subset.var_names))
    matrix = dense(ranking_subset[:, genes].X)
    gene_mean = matrix.mean(axis=0)
    gene_sd = matrix.std(axis=0)
    gene_sd[gene_sd == 0.0] = 1.0
    scaled_cells = (matrix - gene_mean) / gene_sd
    means = []
    scaled_means = []
    for organ in DISPLAY_ORGANS:
        mask = ranking_subset.obs["organ_group"].astype(str).eq(organ).to_numpy()
        means.append(matrix[mask].mean(axis=0))
        scaled_means.append(scaled_cells[mask].mean(axis=0))
    mean_frame = pd.DataFrame(np.vstack(means).T, index=genes, columns=DISPLAY_ORGANS)
    scaled = pd.DataFrame(np.vstack(scaled_means).T, index=genes, columns=DISPLAY_ORGANS)
    mean_frame.to_csv(SCF_PROJECT_ROOT / "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_06/b_S06b/03_PlotData/S6b_All_NaiveB_mean_expression.csv")
    scaled.to_csv(SCF_PROJECT_ROOT / "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_06/b_S06b/03_PlotData/S6b_All_NaiveB_scaled_expression.csv")
    first = min(20, len(scaled))
    second = min(40, len(scaled))
    plot_heatmap(
        scaled,
        panel / "02_Figures" / "FigureS6b_All_NaiveB_Organ_Programs.pdf",
        panel / "02_Figures" / "FigureS6b_All_NaiveB_Organ_Programs.png",
        vmin=-1.0,
        vmax=1.0,
        colorbar_label="Mean expression in group",
        title="All naive B cells",
        annotations=[
            (0, first - 1, "Lymphocyte differentiation"),
            (first, second - 1, "Positive regulation of immune effector process"),
            (second, len(scaled) - 1, "Positive regulation of B cell proliferation"),
        ],
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--canonical", type=Path, default=CANONICAL_ATLAS)
    args = parser.parse_args()
    output = args.output.resolve()
    canonical = args.canonical.resolve()
    if not canonical.is_file():
        raise FileNotFoundError(canonical)

    panels = {
        "S06a": output / "a_S06a",
        "S06b": output / "b_S06b",
    }
    for panel in panels.values():
        for subdir in ("02_Figures", "03_Source_Data", "03_PlotData"):
            (panel / subdir).mkdir(parents=True, exist_ok=True)

    master = sc.read_h5ad(canonical)
    census = validate_inputs(master)
    build_s6a(master, panels["S06a"])
    build_s6b(master, panels["S06b"])
    census.to_csv(SCF_PROJECT_ROOT / "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_06/a_S06a/03_Source_Data/S6_NaiveB_cell_census.tsv", sep="\t", index=False)
    print(f"S6 analysis complete: {int(census['cell_count'].sum())} naive B cells")


if __name__ == "__main__":
    main()
