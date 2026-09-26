"""Plot marker expression across erythroid cell populations for Supplementary Figure 4."""

from __future__ import annotations

import gzip
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from scipy import sparse


GENE_GROUPS = {
    "35_Early_ERY": [
        "TUBA1B", "HMGB2", "MGST3", "UBE2S", "PTTG1",
        "TUBB4B", "CENPF", "NUSAP1", "CD36", "SMC4",
    ],
    "36_Mid_ERY": [
        "XPO7", "MARCH3", "SLC2A1", "RAPGEF2", "SLC25A37",
        "HMBS", "MAN1A1", "HERC1", "SPTA1", "HECTD4",
    ],
    "37_Late_ERY": [
        "UBB", "ALAS2", "NCOA4", "DCAF12", "MKRN1",
        "ADIPOR1", "SLC25A39", "BPGM", "FAM210B", "SNCA",
    ],
}
GROUP_LEVELS = list(GENE_GROUPS)
GENES = [gene for genes in GENE_GROUPS.values() for gene in genes]
def write_source_matrix(data: ad.AnnData, path: Path) -> None:
    matrix = data.X.toarray() if sparse.issparse(data.X) else np.asarray(data.X)
    metadata = data.obs[["MainID", "Main_Organ", "Gestational_Age_Weeks", "Last_cell_type_num"]].copy()
    metadata.insert(0, "cell_id", data.obs_names.astype(str))
    frame = pd.concat(
        [metadata.reset_index(drop=True), pd.DataFrame(matrix, columns=data.var_names)],
        axis=1,
    )
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        frame.to_csv(handle, index=False)


def main() -> int:
    panel_dir = Path(__file__).resolve().parents[1]
    atlas_path = Path(__file__).resolve().parents[4] / "00_Global_Data/Fetal_Immune_Atlas.h5ad"
    tables_dir = panel_dir / "03_Source_Data"
    plots_dir = panel_dir / "02_Figures"
    for directory in (tables_dir, plots_dir):
        directory.mkdir(parents=True, exist_ok=True)
    atlas = ad.read_h5ad(atlas_path, backed="r")
    try:
        selected = atlas.obs["Last_cell_type_num"].astype(str).isin(GROUP_LEVELS)
        data = atlas[selected, GENES].to_memory()
    finally:
        atlas.file.close()
    data.obs["Last_cell_type_num"] = pd.Categorical(
        data.obs["Last_cell_type_num"].astype(str), categories=GROUP_LEVELS, ordered=True
    )
    data.obs["predicted_cell_type"] = pd.Categorical(
        data.obs["Last_cell_type_num"].astype(str), categories=GROUP_LEVELS, ordered=True
    )
    data.uns["predicted_cell_type_colors"] = ["#1F77B4", "#FF7F0E", "#2CA02C"]
    matrix = data.X.tocsr() if sparse.issparse(data.X) else np.asarray(data.X)
    summary_rows: list[dict[str, object]] = []
    for group in GROUP_LEVELS:
        group_mask = np.asarray(data.obs["Last_cell_type_num"] == group)
        values = matrix[group_mask]
        means = np.asarray(values.mean(axis=0)).ravel()
        detected = np.asarray((values > 0).mean(axis=0)).ravel()
        for gene, mean_expression, fraction_expressing in zip(GENES, means, detected):
            summary_rows.append(
                {
                    "cell_type": group,
                    "gene_group": next(key for key, genes in GENE_GROUPS.items() if gene in genes),
                    "gene": gene,
                    "cells": int(group_mask.sum()),
                    "mean_log_normalized_expression": float(mean_expression),
                    "fraction_expressing": float(fraction_expressing),
                }
            )
    pd.DataFrame(summary_rows).to_csv(
        Path(__file__).resolve().parents[4] / "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_04/a_S04a/03_PlotData/S4_expression_summary.csv", index=False
    )
    write_source_matrix(data, tables_dir / "S4_cell_expression_source_data.csv.gz")

    sc.settings.set_figure_params(dpi=120, dpi_save=300, facecolor="white", fontsize=9)
    sc.tl.dendrogram(
        data,
        groupby="predicted_cell_type",
        var_names=GENES,
        use_raw=False,
        cor_method="pearson",
        linkage_method="complete",
        optimal_ordering=False,
        key_added="dendrogram_S4_erythroid",
    )
    axes = sc.pl.heatmap(
        data,
        var_names=GENE_GROUPS,
        groupby="predicted_cell_type",
        use_raw=False,
        standard_scale="var",
        dendrogram="dendrogram_S4_erythroid",
        swap_axes=False,
        cmap="RdBu_r",
        show_gene_labels=True,
        figsize=(7.78, 7.22),
        show=False,
    )
    figure = axes["heatmap_ax"].figure
    figure.set_size_inches(7.78, 7.22)
    pdf_path = plots_dir / "S4_erythropoiesis_heatmap.pdf"
    png_path = plots_dir / "S4_erythropoiesis_heatmap.png"
    figure.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    figure.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(figure)

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
