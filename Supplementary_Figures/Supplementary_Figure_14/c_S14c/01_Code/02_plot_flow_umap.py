"""Draw the spectral-flow UMAP panels from saved event coordinates."""

from __future__ import annotations

from pathlib import Path

SCF_PROJECT_ROOT = Path(__file__).resolve().parents[4]
SCF_GLOBAL = SCF_PROJECT_ROOT / "00_Global_Data"
SCF_SHARED = SCF_GLOBAL / "Shared_Inputs"
from typing import List, Mapping, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D



PBMC_LABEL_OFFSETS: Mapping[int, Tuple[float, float]] = {
    11: (-0.58, 0.62),
    12: (0.58, -0.54),
    14: (-0.62, 0.52),
    15: (0.62, -0.52),
}


CELL_TYPES = [
    "CXCR5- B",
    "CXCR5+ B",
    "Immature NK",
    "Mature NK",
    "Naive CD8+ T",
    "Central memory-like CD8+ T",
    "Effective memory-like CD8+ T",
    "Terminal differentiated CD8+ T",
    "Naive CD4+ T",
    "Central memory-like CD4+ T",
    "Effective memory-like CD4+ T",
    "Terminal differentiated CD4+ T",
    "DNT",
    "LDP",
    "EDP",
    "CD8+ ISP",
    "CD4+ ISP",
    "Ungate",
]

DISPLAY_LABELS = {
    "Naive CD8+ T": "Naïve CD8⁺ T",
    "Central memory-like CD8+ T": "Central memory-like CD8⁺ T",
    "Effective memory-like CD8+ T": "Effective memory-like CD8⁺ T",
    "Terminal differentiated CD8+ T": "Terminal differentiated CD8⁺ T",
    "Naive CD4+ T": "Naïve CD4⁺ T",
    "Central memory-like CD4+ T": "Central memory-like CD4⁺ T",
    "Effective memory-like CD4+ T": "Effective memory-like CD4⁺ T",
    "Terminal differentiated CD4+ T": "Terminal differentiated CD4⁺ T",
    "CD8+ ISP": "CD8⁺ ISP",
    "CD4+ ISP": "CD4⁺ ISP",
}

PALETTE = {
    "CXCR5- B": "#000f99",
    "CXCR5+ B": "#c386ff",
    "Immature NK": "#d66b00",
    "Mature NK": "#cc0047",
    "Naive CD8+ T": "#009999",
    "Central memory-like CD8+ T": "#9900cc",
    "Effective memory-like CD8+ T": "#cc9999",
    "Terminal differentiated CD8+ T": "#ffff00",
    "Naive CD4+ T": "#86c3ff",
    "Central memory-like CD4+ T": "#666666",
    "Effective memory-like CD4+ T": "#ff48ff",
    "Terminal differentiated CD4+ T": "#1900ff",
    "DNT": "#56ac00",
    "LDP": "#336600",
    "EDP": "#ff9933",
    "CD8+ ISP": "#00ccff",
    "CD4+ ISP": "#66ff33",
    "Ungate": "#cccccc",
}

def coordinate_limits(data: pd.DataFrame) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    xmin, xmax = np.quantile(data["UMAP1"], [0.001, 0.999])
    ymin, ymax = np.quantile(data["UMAP2"], [0.001, 0.999])
    xpad = 0.08 * (xmax - xmin)
    ypad = 0.08 * (ymax - ymin)
    return (
        (float(xmin - xpad), float(xmax + xpad)),
        (float(ymin - ypad), float(ymax + ypad)),
    )


def scatter_embedding(
    ax,
    data: pd.DataFrame,
    title: str,
    limits: Tuple[Tuple[float, float], Tuple[float, float]],
    cell_types: Sequence[str],
    palette: Mapping[str, str],
    number_labels: bool = False,
    point_size: float = 1.2,
    label_offsets: Mapping[int, Tuple[float, float]] | None = None,
) -> None:
    for cell_type in reversed(cell_types):
        subset = data.loc[data["cell_type"].eq(cell_type)]
        if subset.empty:
            continue
        ax.scatter(
            subset["UMAP1"],
            subset["UMAP2"],
            s=point_size,
            c=palette[cell_type],
            linewidths=0,
            alpha=0.9 if cell_type != "Ungate" else 0.55,
            rasterized=True,
        )
    if number_labels:
        offsets = label_offsets or {}
        for idx, cell_type in enumerate(cell_types, start=1):
            subset = data.loc[data["cell_type"].eq(cell_type)]
            if subset.empty:
                continue
            x = float(subset["UMAP1"].median())
            y = float(subset["UMAP2"].median())
            dx, dy = offsets.get(idx, (0.0, 0.0))
            annotation_kwargs = {}
            if dx != 0.0 or dy != 0.0:
                annotation_kwargs["arrowprops"] = {
                    "arrowstyle": "-",
                    "color": "#555555",
                    "linewidth": 0.45,
                    "shrinkA": 3,
                    "shrinkB": 3,
                }
            ax.annotate(
                str(idx),
                xy=(x, y),
                xytext=(x + dx, y + dy),
                ha="center",
                va="center",
                fontsize=6.5,
                color="black",
                bbox={
                    "boxstyle": "circle,pad=0.18",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.9,
                },
                zorder=20,
                **annotation_kwargs,
            )
    ax.set_title(title, fontsize=10, pad=4)
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)
        spine.set_color("#222222")


def legend_handles(
    cell_types: Sequence[str],
    palette: Mapping[str, str],
    display_labels: Mapping[str, str],
) -> List[Line2D]:
    return [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markersize=6,
            markerfacecolor=palette[cell_type],
            markeredgecolor="none",
            label=f"{idx}  {display_labels.get(cell_type, cell_type)}",
        )
        for idx, cell_type in enumerate(cell_types, start=1)
    ]


def save_figure(fig, pdf_path: Path, png_path: Path, png_dpi: int = 300) -> None:
    if pdf_path.name.startswith("FigureS14"):
        letter = pdf_path.name[len("FigureS14")]
        target = SCF_PROJECT_ROOT / "Supplementary_Figures/Supplementary_Figure_14" / f"{letter}_S14{letter}" / "02_Figures"
        pdf_path = target / pdf_path.name
        png_path = target / png_path.name
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    fig.savefig(png_path, dpi=png_dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_figures(
    pbmc: pd.DataFrame,
    organ: pd.DataFrame,
    main_dir: Path,
    supp_dir: Path,
    cell_types: Sequence[str],
    palette: Mapping[str, str],
    display_labels: Mapping[str, str],
) -> List[Path]:
    outputs: List[Path] = []
    pbmc_limits = coordinate_limits(pbmc)
    organ_limits = coordinate_limits(organ)

    fig = plt.figure(figsize=(13.2, 5.5), constrained_layout=True)
    grid = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.55])
    ax_pbmc = fig.add_subplot(grid[0, 0])
    ax_organ = fig.add_subplot(grid[0, 1])
    ax_legend = fig.add_subplot(grid[0, 2])
    scatter_embedding(
        ax_pbmc,
        pbmc,
        "PBMC",
        pbmc_limits,
        cell_types,
        palette,
        number_labels=True,
        point_size=1.3,
        label_offsets=PBMC_LABEL_OFFSETS,
    )
    scatter_embedding(
        ax_organ,
        organ,
        "Organ",
        organ_limits,
        cell_types,
        palette,
        number_labels=True,
        point_size=1.1,
    )
    ax_legend.axis("off")
    ax_legend.legend(
        handles=legend_handles(cell_types, palette, display_labels),
        loc="center left",
        frameon=False,
        ncol=2,
        fontsize=8.2,
        handletextpad=0.4,
        columnspacing=1.0,
    )
    main_pdf = main_dir / "02_Figures/Figure6a_Spectral_Flow_UMAP.pdf"
    main_png = main_dir / "02_Figures/Figure6a_Spectral_Flow_UMAP.png"
    save_figure(fig, main_pdf, main_png, png_dpi=400)
    outputs.extend([main_pdf, main_png])

    fig, ax = plt.subplots(figsize=(5.0, 5.0), constrained_layout=True)
    scatter_embedding(
        ax,
        pbmc,
        "PBMC",
        pbmc_limits,
        cell_types,
        palette,
        number_labels=True,
        point_size=1.4,
        label_offsets=PBMC_LABEL_OFFSETS,
    )
    s14a_pdf = supp_dir / "02_Figures/FigureS14a_PBMC_UMAP.pdf"
    s14a_png = supp_dir / "02_Figures/FigureS14a_PBMC_UMAP.png"
    save_figure(fig, s14a_pdf, s14a_png, png_dpi=400)
    outputs.extend([s14a_pdf, s14a_png])

    fig, axes = plt.subplots(1, 4, figsize=(13.0, 3.5), constrained_layout=True)
    scatter_embedding(
        axes[0],
        organ,
        "Organ",
        organ_limits,
        cell_types,
        palette,
        number_labels=True,
        point_size=1.0,
    )
    for ax, organ_name in zip(axes[1:], ["Liver", "Thymus", "Spleen"]):
        scatter_embedding(
            ax,
            organ.loc[organ["organ"].eq(organ_name)],
            organ_name,
            organ_limits,
            cell_types,
            palette,
            point_size=1.2,
        )
    s14b_pdf = supp_dir / "02_Figures/FigureS14b_Organ_UMAP_by_tissue.pdf"
    s14b_png = supp_dir / "02_Figures/FigureS14b_Organ_UMAP_by_tissue.png"
    save_figure(fig, s14b_pdf, s14b_png, png_dpi=400)
    outputs.extend([s14b_pdf, s14b_png])

    pbmc_samples = (
        pbmc[["sample_id", "pcw"]].drop_duplicates().sort_values(["pcw", "sample_id"])
    )
    fig, axes = plt.subplots(4, 5, figsize=(13.5, 10.8), constrained_layout=True)
    for ax, (_, sample_row) in zip(axes.flat, pbmc_samples.iterrows()):
        subset = pbmc.loc[pbmc["sample_id"].eq(sample_row["sample_id"])]
        scatter_embedding(
            ax,
            subset,
            f"{sample_row['pcw']:g} pcw",
            pbmc_limits,
            cell_types,
            palette,
            point_size=1.9,
        )
    for ax in axes.flat[len(pbmc_samples) :]:
        ax.axis("off")
    fig.suptitle("PBMC", fontsize=12)
    s14c_pdf = supp_dir / "02_Figures/FigureS14c_PBMC_UMAP_by_sample.pdf"
    s14c_png = supp_dir / "02_Figures/FigureS14c_PBMC_UMAP_by_sample.png"
    save_figure(fig, s14c_pdf, s14c_png, png_dpi=300)
    outputs.extend([s14c_pdf, s14c_png])

    organ_samples = organ[["sample_id", "organ", "pcw"]].drop_duplicates().copy()
    organ_samples["organ_order"] = organ_samples["organ"].map(
        {"Liver": 0, "Thymus": 1, "Spleen": 2}
    )
    organ_samples = organ_samples.sort_values(["organ_order", "pcw", "sample_id"])
    fig, axes = plt.subplots(3, 5, figsize=(13.5, 8.2), constrained_layout=True)
    for row_idx, organ_name in enumerate(["Liver", "Thymus", "Spleen"]):
        organ_rows = organ_samples.loc[organ_samples["organ"].eq(organ_name)]
        for col_idx, (_, sample_row) in enumerate(organ_rows.iterrows()):
            ax = axes[row_idx, col_idx]
            subset = organ.loc[organ["sample_id"].eq(sample_row["sample_id"])]
            scatter_embedding(
                ax,
                subset,
                f"{sample_row['pcw']:g} pcw\n{sample_row['sample_id']}",
                organ_limits,
                cell_types,
                palette,
                point_size=1.6,
            )
        for col_idx in range(len(organ_rows), 5):
            axes[row_idx, col_idx].axis("off")
        axes[row_idx, 0].set_ylabel(organ_name, fontsize=10, labelpad=8)
    s14d_pdf = supp_dir / "02_Figures/FigureS14d_Organ_UMAP_by_sample.pdf"
    s14d_png = supp_dir / "02_Figures/FigureS14d_Organ_UMAP_by_sample.png"
    save_figure(fig, s14d_pdf, s14d_png, png_dpi=300)
    outputs.extend([s14d_pdf, s14d_png])
    return outputs


def main() -> None:

    source_dir = SCF_SHARED / "Main_Figure_06/a_M06a"
    main_dir = SCF_PROJECT_ROOT / "Main_Figure_06/a_M06a"
    supp_dir = SCF_PROJECT_ROOT / "Supplementary_Figures/Supplementary_Figure_14"
    pbmc_path = source_dir / "M06a_S14_PBMC_event_coordinates.csv.gz"
    organ_path = source_dir / "M06a_S14_Organ_event_coordinates.csv.gz"

    pbmc = pd.read_csv(pbmc_path, compression="gzip")
    organ = pd.read_csv(organ_path, compression="gzip")
    display_labels = {
        key: value.replace("⁺", r"$^{+}$")
        for key, value in DISPLAY_LABELS.items()
    }
    figure_files = render_figures(
        pbmc,
        organ,
        main_dir,
        supp_dir,
        CELL_TYPES,
        PALETTE,
        display_labels,
    )
    print(f"Saved {len(figure_files)} figure files.")


if __name__ == "__main__":
    main()
