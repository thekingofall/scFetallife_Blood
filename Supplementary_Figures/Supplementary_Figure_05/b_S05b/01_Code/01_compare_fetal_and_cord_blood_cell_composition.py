"""Analyze cell composition in Supplementary Figure 5 a-d from fetal and cord blood metadata."""

from __future__ import annotations

import re
from pathlib import Path

SCF_PROJECT_ROOT = Path(__file__).resolve().parents[4]
SCF_GLOBAL = SCF_PROJECT_ROOT / "00_Global_Data"
SCF_SHARED = SCF_GLOBAL / "Shared_Inputs"
from typing import Iterable

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


INPUT_H5AD = SCF_SHARED / "Supplementary_Figures/Supplementary_Figure_05/01_Cell_Metadata/S05_reference_metadata_umap.h5ad"

OUT = Path(__file__).resolve().parents[1]
FIG = OUT / "02_Figures"
SOURCE = SCF_PROJECT_ROOT / "Supplementary_Figures/Supplementary_Figure_05/a_S05a/03_PlotData"

WEEK_THRESHOLD = 26.0

FETAL_COLOR = "#FF8C1A"
CORD_COLOR = "#377EB8"
LATE_COLOR = "#E64B35"
CORD_ENRICH_COLOR = "#4DBBD5"
GREY_COLOR = "#D5D5D5"

plt.rcParams.update(
    {
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }
)


def ensure_directories() -> None:
    for directory in (FIG, SOURCE):
        directory.mkdir(parents=True, exist_ok=True)


def numeric_key(value: str) -> tuple[int, str]:
    match = re.match(r"^(\d+)_", str(value))
    return (int(match.group(1)) if match else 999, str(value))


def sample_week(sample: str) -> float:
    match = re.match(r"^B(\d+(?:\.\d+)?)_", str(sample))
    return float(match.group(1)) if match else np.nan


def stage_for_sample(sample: str) -> str:
    sample = str(sample)
    if sample.startswith("F"):
        return "Cord"
    week = sample_week(sample)
    if np.isnan(week):
        return "Unknown"
    return "Early" if week <= WEEK_THRESHOLD else "Late"


def clean_cell_label(value: str) -> str:
    value = str(value)
    match = re.match(r"^(\d+)_(.*)$", value)
    if not match:
        return value.replace("_", " ")
    number, name = match.groups()
    substitutions = {
        "HSC_MPP": "HSC/MPP",
        "Tem": "Tem-like cells",
        "Th17like_INNATE_T": "Th17-like Innate T",
        "NK T": "NKT",
        "Gamma Delta V1 T": "Vδ1",
        "Gamma Delta V2 T": "Vδ2",
        "GNG4 +CD8aa+T": "GNG4⁺CD8aa⁺ T",
        "Myeloid-CD177": "Myeloid–CD177",
        "CD14+PPBP+ Monocytes": "CD14⁺PPBP⁺ Monocytes",
        "40_Cord": "Cord",
    }
    name = substitutions.get(name, name.replace("_", " "))
    return f"{number}  {name}"


def bh_adjust(p_values: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(p_values), dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)
    if not valid.any():
        return result
    p = values[valid]
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    result[valid] = restored
    return result


def significance_label(p_value: float) -> str:
    if not np.isfinite(p_value):
        return ""
    if p_value < 0.001:
        return "***"
    if p_value < 0.01:
        return "**"
    if p_value < 0.05:
        return "*"
    return ""


def add_umap_arrows(ax: plt.Axes, fontsize: float = 6.0) -> None:
    origin = (0.055, 0.065)
    ax.annotate(
        "",
        xy=(0.16, origin[1]),
        xytext=origin,
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="-|>", lw=0.8, color="black"),
    )
    ax.annotate(
        "",
        xy=(origin[0], 0.18),
        xytext=origin,
        xycoords="axes fraction",
        arrowprops=dict(arrowstyle="-|>", lw=0.8, color="black"),
    )
    ax.text(0.105, 0.012, "UMAP1", transform=ax.transAxes, ha="center", va="bottom", fontsize=fontsize, style="italic")
    ax.text(0.006, 0.12, "UMAP2", transform=ax.transAxes, ha="left", va="center", fontsize=fontsize, rotation=90, style="italic")


def style_umap_axis(ax: plt.Axes, coords: np.ndarray, arrows: bool = False, arrow_fontsize: float = 6.0) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlim(float(coords[:, 0].min()) - 0.5, float(coords[:, 0].max()) + 0.5)
    ax.set_ylim(float(coords[:, 1].min()) - 0.5, float(coords[:, 1].max()) + 0.5)
    ax.set_aspect("equal", adjustable="box")
    if arrows:
        add_umap_arrows(ax, fontsize=arrow_fontsize)


def scatter_colored(
    ax: plt.Axes,
    coords: np.ndarray,
    color_values: np.ndarray,
    size: float,
    alpha: float,
) -> None:
    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=color_values,
        s=size,
        alpha=alpha,
        linewidths=0,
        rasterized=True,
    )


def save_figure(fig: plt.Figure, stem: str, width: float, height: float, dpi: int = 300) -> None:
    fig.set_size_inches(width, height)
    match = re.match(r"S5([a-d])_", stem)
    figure_dir = FIG
    if match:
        letter = match.group(1)
        figure_dir = SCF_PROJECT_ROOT / "Supplementary_Figures/Supplementary_Figure_05" / f"{letter}_S05{letter}" / "02_Figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / f"{stem}.pdf", dpi=dpi)
    fig.savefig(figure_dir / f"{stem}.png", dpi=dpi)
    plt.close(fig)


def build_composition(obs: pd.DataFrame, cell_types: list[str]) -> pd.DataFrame:
    samples = (
        obs[["sample", "week", "stage", "source_label"]]
        .drop_duplicates("sample")
        .assign(sort_week=lambda frame: frame["week"].fillna(999.0))
        .sort_values(["sort_week", "sample"])
        .drop(columns="sort_week")
    )
    index = pd.MultiIndex.from_product(
        [samples["sample"].tolist(), cell_types], names=["sample", "cell_type"]
    )
    counts = obs.groupby(["sample", "harmonized_cell_type"], observed=True).size()
    counts.index = counts.index.set_names(["sample", "cell_type"])
    complete = counts.reindex(index, fill_value=0).rename("cell_count").reset_index()
    complete = complete.merge(samples, on="sample", how="left", validate="many_to_one")
    complete["sample_total"] = complete.groupby("sample")["cell_count"].transform("sum")
    complete["proportion"] = complete["cell_count"] / complete["sample_total"]
    complete["percent"] = complete["proportion"] * 100.0
    complete["cell_type_order"] = complete["cell_type"].map(lambda value: numeric_key(value)[0])
    return complete.sort_values(["week", "sample", "cell_type_order"], na_position="last")


def calculate_enrichment(obs: pd.DataFrame, ordered_cell_types: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    subset = obs.loc[obs["stage"].isin(["Late", "Cord"]), ["stage", "harmonized_cell_type"]].copy()
    table = pd.crosstab(subset["harmonized_cell_type"], subset["stage"], dropna=True)
    for group in ("Late", "Cord"):
        if group not in table.columns:
            table[group] = 0
    table = table[["Late", "Cord"]]
    table = table.loc[table.sum(axis=1) > 0]
    ordered = [cell_type for cell_type in ordered_cell_types if cell_type in table.index]
    table = table.reindex(ordered)
    late_total = int(table["Late"].sum())
    cord_total = int(table["Cord"].sum())
    rows: list[dict[str, object]] = []
    for cell_type, values in table.iterrows():
        late = int(values["Late"])
        cord = int(values["Cord"])
        late_other = late_total - late
        cord_other = cord_total - cord
        fisher_or, p_value = fisher_exact(
            [[late, cord], [late_other, cord_other]], alternative="two-sided"
        )
        log2_or = np.log2(
            ((late + 0.5) * (cord_other + 0.5))
            / ((cord + 0.5) * (late_other + 0.5))
        )
        rows.append(
            {
                "cell_type": cell_type,
                "late_cells": late,
                "cord_cells": cord,
                "late_other_cells": late_other,
                "cord_other_cells": cord_other,
                "late_total": late_total,
                "cord_total": cord_total,
                "fisher_odds_ratio": fisher_or,
                "log2_odds_ratio_pseudocount_0.5": log2_or,
                "p_value_raw": p_value,
            }
        )
    result = pd.DataFrame(rows)
    result["q_value_BH"] = bh_adjust(result["p_value_raw"])
    result["significance_raw_p"] = result["p_value_raw"].map(significance_label)
    result["significance_BH"] = result["q_value_BH"].map(significance_label)
    table.index.name = "cell_type"
    return table, result


def a_legend_handles(categories: list[str], palette: dict[str, str], markersize: float) -> list[Line2D]:
    return [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=palette[category],
            markeredgecolor="none",
            label=clean_cell_label(category),
            markersize=markersize,
        )
        for category in categories
    ]


def plot_panel_a(
    coords: np.ndarray,
    display_types: pd.Series,
    categories: list[str],
    palette: dict[str, str],
) -> None:
    fig = plt.figure(figsize=(10.3, 3.5))
    ax = fig.add_axes([0.035, 0.08, 0.42, 0.86])
    colors = display_types.map(palette).to_numpy()
    scatter_colored(ax, coords, colors, size=0.10, alpha=0.80)
    style_umap_axis(ax, coords, arrows=True, arrow_fontsize=6.5)
    handles = a_legend_handles(categories, palette, markersize=4.7)
    fig.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(0.46, 0.48),
        frameon=False,
        ncol=3,
        columnspacing=1.5,
        handletextpad=0.35,
        fontsize=6.9,
    )
    fig.text(0.705, 0.95, "Cell type", ha="center", va="top", fontsize=9.5, fontweight="bold")
    save_figure(fig, "S5a_integrated_UMAP_40Cord", 10.3, 3.5)


def plot_panel_b(coords: np.ndarray, source_labels: pd.Series) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.0))
    groups = [
        ("Cord Blood (GSE157007)", "Cord Blood (GSE157007)", CORD_COLOR),
        ("Fetal Blood in this study", "Fetal Blood in this study", FETAL_COLOR),
    ]
    for index, (group, title, color) in enumerate(groups):
        ax = axes[index]
        mask = source_labels.eq(group).to_numpy()
        ax.scatter(
            coords[:, 0], coords[:, 1], c=GREY_COLOR, s=0.075, alpha=0.50,
            linewidths=0, rasterized=True,
        )
        ax.scatter(
            coords[mask, 0], coords[mask, 1], c=color, s=0.12, alpha=0.85,
            linewidths=0, rasterized=True,
        )
        style_umap_axis(ax, coords, arrows=index == 0, arrow_fontsize=5.5)
        ax.set_title(title, fontsize=9.2, fontweight="bold", pad=2)
    fig.subplots_adjust(left=0.03, right=0.99, bottom=0.05, top=0.90, wspace=0.14)
    save_figure(fig, "S5b_dataset_small_multiples", 7.4, 3.0)


def composition_plot_data(
    composition: pd.DataFrame, cell_types: list[str]
) -> tuple[list[str], pd.DataFrame]:
    samples = (
        composition[["sample", "week", "stage"]]
        .drop_duplicates("sample")
        .assign(sort_week=lambda frame: frame["week"].fillna(999.0))
        .sort_values(["sort_week", "sample"])["sample"]
        .tolist()
    )
    wide = composition.pivot(index="sample", columns="cell_type", values="proportion").fillna(0.0)
    wide = wide.reindex(index=samples, columns=cell_types, fill_value=0.0)
    return samples, wide


def draw_composition(
    ax: plt.Axes,
    wide: pd.DataFrame,
    cell_types: list[str],
    palette: dict[str, str],
    xtick_fontsize: float,
    ytick_fontsize: float,
    ylabel_fontsize: float,
) -> None:
    x = np.arange(len(wide))
    bottom = np.zeros(len(wide), dtype=float)
    for cell_type in cell_types:
        values = wide[cell_type].to_numpy(dtype=float)
        ax.bar(x, values, bottom=bottom, color=palette[cell_type], width=0.82, linewidth=0)
        bottom += values
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Composition (%)", fontsize=ylabel_fontsize)
    ticks = np.linspace(0.0, 1.0, 6)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{int(value * 100)}%" for value in ticks], fontsize=ytick_fontsize)
    ax.set_xticks(x)
    ax.set_xticklabels(wide.index.tolist(), rotation=67, ha="right", fontsize=xtick_fontsize)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.margins(x=0.005)


def add_composition_group_headers(ax: plt.Axes, wide: pd.DataFrame, fontsize: float) -> None:
    fetal_count = sum(str(sample).startswith("B") for sample in wide.index)
    cord_count = len(wide) - fetal_count
    if fetal_count:
        left, right = -0.40, fetal_count - 0.60
        ax.plot([left, right], [1.03, 1.03], color="black", lw=1.0, clip_on=False)
        ax.text((left + right) / 2, 1.06, "Fetal Blood in this study", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")
    if cord_count:
        left, right = fetal_count - 0.40, len(wide) - 0.60
        ax.plot([left, right], [1.03, 1.03], color="black", lw=1.0, clip_on=False)
        ax.text((left + right) / 2, 1.06, "Cord Blood\n(GSE157007)", ha="center", va="bottom", fontsize=fontsize, fontweight="bold")


def plot_panel_c(
    composition: pd.DataFrame,
    cell_types: list[str],
    palette: dict[str, str],
) -> None:
    _, wide = composition_plot_data(composition, cell_types)
    fig = plt.figure(figsize=(15.0, 6.0))
    ax = fig.add_axes([0.055, 0.25, 0.69, 0.64])
    draw_composition(ax, wide, cell_types, palette, 7.2, 8.2, 9.5)
    add_composition_group_headers(ax, wide, fontsize=8.5)
    handles = [Patch(facecolor=palette[cell_type], edgecolor="none", label=clean_cell_label(cell_type)) for cell_type in cell_types]
    fig.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(0.76, 0.54),
        frameon=False,
        ncol=2,
        columnspacing=1.0,
        handletextpad=0.4,
        fontsize=6.8,
    )
    fig.text(0.875, 0.92, "Cell type", ha="center", va="top", fontsize=9.5, fontweight="bold")
    save_figure(fig, "S5c_sample_composition", 15.0, 6.0)


def draw_enrichment(
    ax: plt.Axes,
    enrichment: pd.DataFrame,
    tick_fontsize: float,
    label_fontsize: float,
    title_fontsize: float,
    star_fontsize: float,
) -> None:
    values = enrichment["log2_odds_ratio_pseudocount_0.5"].to_numpy(dtype=float)
    x = np.arange(len(enrichment))
    colors = np.where(values > 0, LATE_COLOR, CORD_ENRICH_COLOR)
    ax.bar(x, values, color=colors, width=0.82, linewidth=0)
    ax.axhline(0.0, color="#777777", linewidth=0.7)
    ax.set_ylim(-3.0, 10.0)
    ax.set_ylabel("Enrichment (log2OR)", fontsize=label_fontsize)
    ax.set_title("Enrichment: Late vs Cord", fontsize=title_fontsize, pad=5)
    ax.set_xticks(x)
    ax.set_xticklabels(
        enrichment["cell_type"].tolist(), rotation=90, ha="center", fontsize=tick_fontsize
    )
    ax.tick_params(axis="y", labelsize=tick_fontsize)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for index, row in enrichment.reset_index(drop=True).iterrows():
        star = row["significance_raw_p"]
        if not star:
            continue
        value = float(row["log2_odds_ratio_pseudocount_0.5"])
        y = min(value + 0.18, 9.65) if value >= 0 else max(value - 0.18, -2.65)
        ax.text(
            index,
            y,
            star,
            ha="center",
            va="bottom" if value >= 0 else "top",
            fontsize=star_fontsize,
        )
    ax.margins(x=0.008)


def enrichment_legend_handles() -> list[Patch]:
    return [
        Patch(facecolor=LATE_COLOR, edgecolor="none", label="Fetal Blood (>26pcw) in this study"),
        Patch(facecolor=CORD_ENRICH_COLOR, edgecolor="none", label="Cord Blood (GSE157007)"),
    ]


def plot_panel_d(enrichment: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(14.0, 5.3))
    ax = fig.add_axes([0.055, 0.33, 0.70, 0.58])
    draw_enrichment(ax, enrichment, 7.1, 9.8, 10.5, 7.8)
    fig.legend(
        handles=enrichment_legend_handles(),
        title="Group",
        loc="center left",
        bbox_to_anchor=(0.78, 0.64),
        frameon=False,
        fontsize=8.2,
        title_fontsize=8.5,
    )
    save_figure(fig, "S5d_late_vs_cord_enrichment_rawP", 14.0, 5.3)


def plot_composite(
    coords: np.ndarray,
    display_types: pd.Series,
    source_labels: pd.Series,
    a_categories: list[str],
    a_palette: dict[str, str],
    composition: pd.DataFrame,
    c_cell_types: list[str],
    c_palette: dict[str, str],
    enrichment: pd.DataFrame,
) -> None:
    fig = plt.figure(figsize=(8.5, 11.0))
    outer = fig.add_gridspec(
        4,
        1,
        height_ratios=[2.30, 2.15, 2.85, 2.55],
        left=0.055,
        right=0.985,
        bottom=0.085,
        top=0.98,
        hspace=0.30,
    )


    grid_a = outer[0].subgridspec(1, 2, width_ratios=[1.0, 1.50], wspace=0.03)
    ax_a = fig.add_subplot(grid_a[0])
    scatter_colored(ax_a, coords, display_types.map(a_palette).to_numpy(), size=0.06, alpha=0.80)
    style_umap_axis(ax_a, coords, arrows=True, arrow_fontsize=4.1)
    legend_a = fig.add_subplot(grid_a[1])
    legend_a.axis("off")
    legend_a.legend(
        handles=a_legend_handles(a_categories, a_palette, 2.8),
        loc="center",
        frameon=False,
        ncol=3,
        columnspacing=0.65,
        handletextpad=0.20,
        fontsize=4.0,
        title="Cell type",
        title_fontsize=5.5,
    )


    grid_b = outer[1].subgridspec(1, 2, wspace=0.12)
    b_groups = [
        ("Cord Blood (GSE157007)", CORD_COLOR),
        ("Fetal Blood in this study", FETAL_COLOR),
    ]
    for index, (group, color) in enumerate(b_groups):
        ax = fig.add_subplot(grid_b[index])
        mask = source_labels.eq(group).to_numpy()
        ax.scatter(coords[:, 0], coords[:, 1], c=GREY_COLOR, s=0.045, alpha=0.45, linewidths=0, rasterized=True)
        ax.scatter(coords[mask, 0], coords[mask, 1], c=color, s=0.065, alpha=0.85, linewidths=0, rasterized=True)
        style_umap_axis(ax, coords, arrows=index == 0, arrow_fontsize=3.7)
        ax.set_title(group, fontsize=5.3, fontweight="bold", pad=1)


    grid_c = outer[2].subgridspec(1, 2, width_ratios=[2.90, 1.0], wspace=0.02)
    ax_c = fig.add_subplot(grid_c[0])
    _, wide = composition_plot_data(composition, c_cell_types)
    draw_composition(ax_c, wide, c_cell_types, c_palette, 4.0, 4.5, 5.3)
    add_composition_group_headers(ax_c, wide, fontsize=4.7)
    legend_c = fig.add_subplot(grid_c[1])
    legend_c.axis("off")
    legend_c.legend(
        handles=[Patch(facecolor=c_palette[cell_type], edgecolor="none", label=clean_cell_label(cell_type)) for cell_type in c_cell_types],
        loc="center",
        frameon=False,
        ncol=2,
        columnspacing=0.55,
        handletextpad=0.25,
        fontsize=3.2,
        title="Cell type",
        title_fontsize=4.8,
    )


    grid_d = outer[3].subgridspec(1, 2, width_ratios=[3.30, 1.0], wspace=0.04)
    ax_d = fig.add_subplot(grid_d[0])
    draw_enrichment(ax_d, enrichment, 3.8, 5.0, 5.8, 4.2)
    legend_d = fig.add_subplot(grid_d[1])
    legend_d.axis("off")
    legend_d.legend(
        handles=enrichment_legend_handles(),
        title="Group",
        loc="center",
        frameon=False,
        fontsize=4.2,
        title_fontsize=4.8,
    )

    for label, y in (("a", 0.974), ("b", 0.738), ("c", 0.505), ("d", 0.245)):
        fig.text(0.012, y, label, fontsize=13, fontweight="bold", ha="left", va="top")

    save_figure(fig, "Supplementary_Figure5_a-d", 8.5, 11.0)


def main() -> int:
    ensure_directories()
    backed = ad.read_h5ad(INPUT_H5AD, backed="r")
    try:
        obs_all = backed.obs.copy()
        coords_no = np.asarray(backed.obsm["X_umap"][:], dtype=np.float32)
        original_colors = [str(v) for v in backed.uns["Last_cell_type_num_colors"]]
    finally:
        backed.file.close()
    if coords_no.shape != (len(obs_all), 2):
        raise ValueError("Unexpected fetal and cord blood atlas dimensions")
    obs_all = obs_all.copy()
    obs_all["sample"] = obs_all["MainID"].astype(str)
    obs_all["week"] = obs_all["sample"].map(sample_week)
    obs_all["stage"] = obs_all["sample"].map(stage_for_sample)
    obs_all["source_label"] = np.where(
        obs_all["stage"].eq("Cord"),
        "Cord Blood (GSE157007)",
        "Fetal Blood in this study",
    )
    obs_all["display_cell_type_S5a"] = np.where(
        obs_all["stage"].eq("Cord"),
        "40_Cord",
        obs_all["Last_cell_type_num"].astype(str),
    )
    obs_all["harmonized_cell_type"] = np.where(
        obs_all["stage"].eq("Cord"),
        obs_all["predicted_cell_type"].astype(str),
        obs_all["Last_cell_type_num"].astype(str),
    )

    obs_no = obs_all
    a_categories = sorted(obs_no["display_cell_type_S5a"].unique().tolist(), key=numeric_key)
    c_cell_types = sorted(obs_no["harmonized_cell_type"].unique().tolist(), key=numeric_key)
    submitted_categories = (
        obs_all["Last_cell_type_num"].cat.categories.astype(str).tolist()
        if isinstance(obs_all["Last_cell_type_num"].dtype, pd.CategoricalDtype)
        else sorted(obs_all["Last_cell_type_num"].unique().astype(str).tolist(), key=numeric_key)
    )
    submitted_categories = [category for category in submitted_categories if category in a_categories]
    submitted_categories = sorted(submitted_categories, key=numeric_key)
    if len(submitted_categories) != len(original_colors):
        raise RuntimeError(
            f"Cell-type category and color counts differ: {len(submitted_categories)} vs {len(original_colors)}"
        )
    a_palette = dict(zip(submitted_categories, original_colors))
    if set(a_categories) != set(a_palette):
        raise RuntimeError("S5a palette is missing cell-type colors")

    cmap_b = plt.get_cmap("tab20b", 20)
    cmap_c = plt.get_cmap("tab20c", 20)
    c_colors = [matplotlib.colors.to_hex(cmap_b(index)) for index in range(20)]
    c_colors.extend(
        matplotlib.colors.to_hex(cmap_c(index)) for index in range(len(c_cell_types) - 20)
    )
    c_palette = dict(zip(c_cell_types, c_colors))

    composition = build_composition(obs_no, c_cell_types)
    contingency, enrichment = calculate_enrichment(obs_no, c_cell_types)
    composition.to_csv(SCF_PROJECT_ROOT / "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/a_S05a/03_PlotData/S5c_composition_complete.csv", index=False)
    contingency.to_csv(SCF_PROJECT_ROOT / "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/a_S05a/03_PlotData/S5d_contingency.csv")
    enrichment.to_csv(SCF_PROJECT_ROOT / "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/a_S05a/03_PlotData/S5d_fisher_rawP_BHsecondary.csv", index=False)
    cell_metadata = obs_no[
        [
            "sample",
            "week",
            "stage",
            "source_label",
            "display_cell_type_S5a",
            "harmonized_cell_type",
            "Last_cell_type_num",
            "predicted_cell_type",
        ]
    ].copy()
    cell_metadata.insert(0, "cell_id", cell_metadata.index.astype(str))
    cell_metadata["UMAP1"] = coords_no[:, 0]
    cell_metadata["UMAP2"] = coords_no[:, 1]
    cell_metadata.to_csv(
        SCF_PROJECT_ROOT / "00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/a_S05a/03_PlotData/S5ab_cell_metadata_coordinates.csv.gz",
        index=False,
        compression="gzip",
    )

    plot_panel_a(coords_no, obs_no["display_cell_type_S5a"], a_categories, a_palette)
    plot_panel_b(coords_no, obs_no["source_label"])
    plot_panel_c(composition, c_cell_types, c_palette)
    plot_panel_d(enrichment)
    plot_composite(coords_no, obs_no["display_cell_type_S5a"], obs_no["source_label"],
                   a_categories, a_palette, composition, c_cell_types, c_palette, enrichment)
    print(f"S5 complete: {len(obs_no)} cells; {obs_no['sample'].nunique()} samples")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
