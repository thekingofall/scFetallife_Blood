"""Calculate spectral-flow UMAP coordinates from gated flow-cytometry events.

FlowJo 10.8.1 workspaces specify the gates applied by FlowKit.
The DownsampleV3 algebraic score selects 5,000 PBMC or
9,000 organ parent-population events per sample (bounded by the available parent
count); and UMAP uses the parameters reported in the manuscript (Euclidean,
15 neighbours, min_dist 0.5, two components).
"""

from __future__ import annotations

import argparse
import gc
import math
import os
import re
import xml.etree.ElementTree as ET
import flowkit as fk
import flowutils
from flowkit._models.transforms._transforms import LinearTransform
from flowkit._utils import wsp_utils as flowkit_wsp_utils
import numpy as np
import pandas as pd
from scipy.linalg import orthogonal_procrustes
import umap

from pathlib import Path

SCF_PROJECT_ROOT = Path(__file__).resolve().parents[4]
SCF_GLOBAL = SCF_PROJECT_ROOT / "00_Global_Data"
SCF_SHARED = SCF_GLOBAL / "Shared_Inputs"
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple




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


GATE_NAMES = {
    "CXCR5- B": ("CXCR5-B",),
    "CXCR5+ B": ("CXCR5+B",),
    "Immature NK": ("Immature NK",),
    "Mature NK": ("Mature NK",),
    "Naive CD8+ T": ("Naive CD8+T",),
    "Central memory-like CD8+ T": ("Central memory-like CD8+T",),
    "Effective memory-like CD8+ T": ("Effective memory-like CD8+T",),
    "Terminal differentiated CD8+ T": ("Terminal differentiated CD8+T",),
    "Naive CD4+ T": ("Naive CD4+T",),
    "Central memory-like CD4+ T": ("Central memory-like CD4+T",),
    "Effective memory-like CD4+ T": ("Effective memory-like CD4+T",),
    "Terminal differentiated CD4+ T": ("Terminal differentiated CD4+T",),
    "DNT": ("CD3+CD4-CD8-", "CD3+DNT"),
    "LDP": ("LDP",),
    "EDP": ("EDP",),
    "CD8+ ISP": ("CD8+ISP",),
    "CD4+ ISP": ("CD4+ISP",),
}


ORGAN_TARGETS = {
    "CXCR5- B": (-2.1, -0.3),
    "CXCR5+ B": (-1.9, -1.5),
    "Immature NK": (-0.8, 2.1),
    "Mature NK": (-0.6, 1.3),
    "Naive CD8+ T": (1.7, 0.4),
    "Naive CD4+ T": (2.0, -0.9),
    "EDP": (0.3, -1.4),
}



def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]




def normalise_flow_id(value: str) -> str:
    match = re.fullmatch(r"([A-Za-z]+)0*([0-9]+)", str(value).strip())
    if not match:
        return str(value).strip().upper()
    return f"{match.group(1).upper()}{int(match.group(2))}"


def flow_id_from_filename(filename: str) -> str:
    stem = Path(filename).stem
    match = re.search(r"[-_]([A-Za-z]+0*[0-9]+)(?:_Unmixed)?$", stem, flags=re.I)
    if not match:
        raise ValueError(f"Cannot parse Flow ID from {filename}")
    return normalise_flow_id(match.group(1))


def read_sample_metadata(path: Path) -> pd.DataFrame:
    table = pd.read_excel(path, sheet_name="Sheet1", header=1)
    table.columns = [str(c).strip() for c in table.columns]
    table["Fetal ID"] = table["Fetal ID"].ffill()
    table["Flow ID normalized"] = table["Flow ID"].map(normalise_flow_id)
    table["Sample ID"] = table["Sample ID"].astype(str).str.strip()
    table["Fetal organ"] = table["Fetal organ"].astype(str).str.strip()
    table["pcw"] = pd.to_numeric(table["pcw"])

    return table


def parse_downsample_specs(wsp_path: Path) -> Dict[str, Dict[str, float]]:
    specs: Dict[str, Dict[str, float]] = {}
    for _, element in ET.iterparse(str(wsp_path), events=("end",)):
        if local_name(element.tag) != "Sample":
            continue
        filename = None
        formula = None
        gate_min = None
        gate_max = None
        target = None
        for node in element.iter():
            lname = local_name(node.tag)
            if lname == "Keyword" and node.attrib.get("name") == "$FIL":
                filename = node.attrib.get("value")
            elif lname == "DerivedParameter" and node.attrib.get("name") == "DownsampleDP":
                formula = node.attrib.get("formula")
            elif lname == "DownSample" and node.attrib.get("numEvents"):
                target = int(node.attrib["numEvents"])
            elif lname == "dimension":
                children = list(node)
                if children and children[0].attrib.get(
                    "{http://www.isac-net.org/std/Gating-ML/v2.0/datatypes}name"
                ) == "DownsampleDP":
                    gate_min = float(
                        node.attrib.get(
                            "{http://www.isac-net.org/std/Gating-ML/v2.0/gating}min",
                            node.attrib.get("min", "nan"),
                        )
                    )
                    gate_max = float(
                        node.attrib.get(
                            "{http://www.isac-net.org/std/Gating-ML/v2.0/gating}max",
                            node.attrib.get("max", "nan"),
                        )
                    )
        if filename and formula:
            mod_match = re.search(r"#%([0-9]+)", formula)
            total_match = re.search(r"/\(([0-9]+)\s*\+\s*1\)", formula)
            if not mod_match or not total_match:
                raise ValueError(f"Unrecognised DownsampleDP formula for {filename}: {formula}")
            specs[Path(filename).name] = {
                "modulus": int(mod_match.group(1)),
                "formula_total": int(total_match.group(1)),
                "gate_min": gate_min,
                "gate_max": gate_max,
                "target": target or 5000,
                "formula": formula,
            }
        element.clear()
    return specs



def gate_id_by_name(gs, candidates: Sequence[str]) -> Tuple[str, Tuple[str, ...]]:
    matches = []
    candidate_set = {c.strip() for c in candidates}
    for gate_name, gate_path in gs.get_gate_ids():
        if gate_name.strip() in candidate_set:
            matches.append((gate_name, gate_path))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one gate for {candidates}, found {matches}")
    return matches[0]


def membership(results, gs, candidates: Sequence[str]) -> np.ndarray:
    gate_name, gate_path = gate_id_by_name(gs, candidates)
    return results.get_gate_membership(gate_name, gate_path=gate_path)


def load_workspace_with_flowjo_linear_minrange_compatibility(
    wsp_path: Path, fcs_paths: Sequence[Path], enable: bool
) -> object:
    """Apply FlowJo linear transforms to the workspace gates and events.

    Convert the lower endpoint m to A=-m for the GatingML formula (x+A)/(T+A).
    """

    if not enable:
        workspace = fk.Workspace(str(wsp_path), fcs_samples=[str(path) for path in fcs_paths])
        return workspace

    original_parser = flowkit_wsp_utils._parse_wsp_transforms

    def corrected_parser(transforms_element, *parser_args, **parser_kwargs):
        transform_lut = original_parser(
            transforms_element, *parser_args, **parser_kwargs
        )
        for channel, transform in transform_lut.items():
            if isinstance(transform, LinearTransform) and float(transform.param_a) < 0:
                old_a = float(transform.param_a)
                transform.param_a = -old_a
        return transform_lut

    flowkit_wsp_utils._parse_wsp_transforms = corrected_parser
    try:
        workspace = fk.Workspace(
            str(wsp_path), fcs_samples=[str(path) for path in fcs_paths]
        )
    finally:
        flowkit_wsp_utils._parse_wsp_transforms = original_parser
    return workspace



def systematic_downsample(
    parent: np.ndarray, spec: Mapping[str, float]
) -> Tuple[np.ndarray, np.ndarray, int]:
    total = len(parent)
    event_number = np.arange(total, dtype=np.float64)
    score = (
        np.mod(event_number, int(spec["modulus"]))
        - 0.5
        + event_number / (float(spec["formula_total"]) + 1.0)
    )
    parent_idx = np.flatnonzero(parent)
    if math.isfinite(float(spec["gate_min"])) and math.isfinite(float(spec["gate_max"])):
        in_saved_gate = (
            (score[parent_idx] >= float(spec["gate_min"]))
            & (score[parent_idx] <= float(spec["gate_max"]))
        )
        saved_gate_count = int(in_saved_gate.sum())
    else:
        saved_gate_count = -1
    target = min(int(spec["target"]), len(parent_idx))
    order = np.argsort(score[parent_idx], kind="mergesort")
    selected_idx = parent_idx[order[:target]]
    return selected_idx, score[selected_idx], saved_gate_count


def compensate_and_transform(sample, gs, selected_idx: np.ndarray) -> Tuple[np.ndarray, List[str]]:
    raw = sample.get_events(source="raw")[selected_idx, :].copy()
    matrix = gs.comp_matrices["Acquisition-defined"]
    detector_indices = [sample.get_channel_index(d) for d in matrix.detectors]
    comp = flowutils.compensate.compensate(raw, matrix.matrix, detector_indices)

    marker_indices = [
        i
        for i, marker in enumerate(sample.pns_labels)
        if str(marker).strip() and str(marker).strip() != "L-D"
    ]
    markers = [str(sample.pns_labels[i]).strip() for i in marker_indices]
    columns = []
    for idx in marker_indices:
        detector = sample.pnn_labels[idx]
        transform_name = f"Comp-{detector}"
        if transform_name not in gs.transformations:
            raise KeyError(f"Missing FlowJo transform {transform_name} for {sample.id}")
        transformed = gs.transformations[transform_name].apply(comp[:, [idx]])
        columns.append(np.asarray(transformed).reshape(-1))
    x = np.column_stack(columns).astype(np.float32, copy=False)
    if not np.isfinite(x).all():
        raise ValueError(f"Non-finite transformed values in {sample.id}")
    return x, markers




def process_compartment(
    compartment: str,
    wsp_path: Path,
    fcs_dir: Path,
    metadata: pd.DataFrame,
    cache_dir: Path,
) -> List[Dict]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    fcs_paths = sorted(fcs_dir.glob("*.fcs"))
    fcs_by_flow: Dict[str, Path] = {}
    for fcs_path in fcs_paths:
        try:
            flow_id = flow_id_from_filename(fcs_path.name)
        except ValueError:
            continue
        if flow_id in fcs_by_flow:
            raise RuntimeError(f"Duplicate FCS Flow ID {flow_id} in {fcs_dir}")
        fcs_by_flow[flow_id] = fcs_path
    workspace = (
        load_workspace_with_flowjo_linear_minrange_compatibility(
            wsp_path, fcs_paths, enable=compartment == "Organ"
        )
    )
    wsp_ids = workspace.get_sample_ids()
    id_by_flow = {}
    for sample_id in wsp_ids:
        try:
            id_by_flow[flow_id_from_filename(sample_id)] = sample_id
        except ValueError:
            continue
    specs = parse_downsample_specs(wsp_path)
    specs_by_flow = {}
    for filename, spec in specs.items():
        try:
            flow_id = flow_id_from_filename(filename)
        except ValueError:
            continue
        specs_by_flow[flow_id] = (filename, spec)

    if compartment == "PBMC":
        target_meta = metadata.loc[metadata["Fetal organ"].eq("PBMC")].copy()
        parent_candidates = ("Live",)
    else:
        target_meta = metadata.loc[metadata["Fetal organ"].isin(["Liver", "Thymus", "Spleen"])].copy()
        parent_candidates = ("CD45+ cell",)

    sample_records: List[Dict] = []
    target_meta = target_meta.sort_values(["pcw", "Fetal organ", "Sample ID"])
    for row in target_meta.to_dict("records"):
        flow_id = row["Flow ID normalized"]
        if flow_id not in id_by_flow:
            raise KeyError(f"No workspace sample for Flow ID {flow_id}")
        if flow_id not in fcs_by_flow:
            raise KeyError(f"No FCS file for Flow ID {flow_id}")
        sample_id = id_by_flow[flow_id]
        fcs_path = fcs_by_flow[flow_id]
        cache_path = cache_dir / f"{row['Sample ID']}.npz"
        print(f"Gate {compartment} {row['Sample ID']} ({sample_id})", flush=True)
        sample = workspace.get_sample(sample_id)
        gs = workspace.get_gating_strategy(sample_id)
        downsample_ids = [
            (name, path)
            for name, path in gs.get_gate_ids()
            if name.strip() == "DownsampleDP"
        ]
        if len(downsample_ids) > 1:
            raise RuntimeError(f"Expected at most one DownsampleDP gate in {sample_id}")
        if downsample_ids:
            gs.remove_gate(downsample_ids[0][0], gate_path=downsample_ids[0][1])
        results = gs.gate_sample(sample)

        parent = membership(results, gs, parent_candidates)
        parent_count = int(parent.sum())
        if parent_count == 0:
            raise ValueError(f"Empty parent population in {sample_id}")
        member_matrix = []
        for cell_type in CELL_TYPES[:-1]:
            gate_membership = membership(results, gs, GATE_NAMES[cell_type])
            member_matrix.append(gate_membership)
        member_matrix_np = np.vstack(member_matrix)

        if flow_id in specs_by_flow:
            spec_filename, spec = specs_by_flow[flow_id]
        else:
            saved_targets = [int(saved_spec[1]["target"]) for saved_spec in specs_by_flow.values()]
            if not saved_targets:
                raise RuntimeError(f"No saved DownsampleV3 target is available for {compartment}")
            target = int(pd.Series(saved_targets).mode().iloc[0])
            modulus = max(1, int(round(sample.event_count / target)))
            spec_filename = sample_id
            spec = {
                "modulus": modulus,
                "formula_total": int(sample.event_count),
                "gate_min": float("nan"),
                "gate_max": float("nan"),
                "target": target,
                "formula": (
                    f"(Event #%{modulus}) - 0.5 + "
                    f"(Event #/({sample.event_count} + 1))"
                ),
            }
        if int(spec["formula_total"]) != sample.event_count:
            raise AssertionError(
                f"Formula event total mismatch for {sample_id}: "
                f"{spec['formula_total']} != {sample.event_count}"
            )
        selected_idx, selected_score, saved_gate_count = systematic_downsample(parent, spec)
        selected_members = member_matrix_np[:, selected_idx]
        membership_sum = selected_members.sum(axis=0)
        label_code = np.full(len(selected_idx), len(CELL_TYPES) - 1, dtype=np.int16)
        for idx in range(len(CELL_TYPES) - 1):
            label_code[selected_members[idx]] = idx

        x, marker_names = compensate_and_transform(sample, gs, selected_idx)
        temporary_cache_path = cache_path.with_suffix(".tmp.npz")
        np.savez_compressed(
            temporary_cache_path,
            X=x,
            marker_names=np.asarray(marker_names, dtype="U64"),
            label_code=label_code,
            event_index=selected_idx.astype(np.int64),
            downsample_score=selected_score.astype(np.float64),
        )
        os.replace(temporary_cache_path, cache_path)

        sample_record = {
            "compartment": compartment,
            "sample_id": row["Sample ID"],
            "flow_id": flow_id,
            "fetal_id": row["Fetal ID"],
            "organ": row["Fetal organ"],
            "pcw": float(row["pcw"]),
            "workspace_sample_id": sample_id,
            "raw_event_count": int(sample.event_count),
            "parent_population": parent_candidates[0],
            "parent_count_flowkit": parent_count,
            "downsample_target": int(spec["target"]),
            "selected_event_count": int(len(selected_idx)),
            "downsample_modulus": int(spec["modulus"]),
            "downsample_formula": spec["formula"],
            "downsample_gate_min": (
                float(spec["gate_min"])
                if math.isfinite(float(spec["gate_min"]))
                else None
            ),
            "downsample_gate_max": (
                float(spec["gate_max"])
                if math.isfinite(float(spec["gate_max"]))
                else None
            ),
            "marker_count": len(marker_names),
            "marker_names": marker_names,
            "cache_path": str(cache_path),
        }
        sample_records.append(sample_record)
        gs.clear_cache()
        del results, gs, sample, x, member_matrix_np
        gc.collect()
    return sample_records


def load_cache_set(
    records: Sequence[Mapping]
) -> Tuple[np.ndarray, pd.DataFrame, List[str]]:
    xs = []
    metadata_rows = []
    marker_names_ref = None
    for record in records:
        with np.load(record["cache_path"], allow_pickle=False) as data:
            markers = [str(x) for x in data["marker_names"]]
            x = data["X"]
            event_index = data["event_index"]
            downsample_score = data["downsample_score"]
            label_code = data["label_code"]
        if marker_names_ref is None:
            marker_names_ref = markers
        elif markers != marker_names_ref:
            raise AssertionError(f"Marker mismatch in {record['cache_path']}")
        xs.append(x)
        for i in range(len(x)):
            metadata_rows.append(
                {
                    "sample_id": record["sample_id"],
                    "flow_id": record["flow_id"],
                    "fetal_id": record["fetal_id"],
                    "organ": record["organ"],
                    "pcw": record["pcw"],
                    "event_index": int(event_index[i]),
                    "downsample_score": float(downsample_score[i]),
                    "cell_type": CELL_TYPES[int(label_code[i])],
                }
            )
    return np.vstack(xs), pd.DataFrame(metadata_rows), marker_names_ref or []


def fit_umap(x: np.ndarray, seed: int) -> np.ndarray:
    model = umap.UMAP(
        n_neighbors=15,
        min_dist=0.5,
        n_components=2,
        metric="euclidean",
        random_state=seed,
        transform_seed=seed,
        n_jobs=1,
        low_memory=True,
        verbose=True,
    )
    return model.fit_transform(x).astype(np.float32)


def align_similarity(source: np.ndarray, target: np.ndarray) -> Tuple[np.ndarray, Dict]:
    source_mean = source.mean(axis=0)
    target_mean = target.mean(axis=0)
    source_centered = source - source_mean
    target_centered = target - target_mean
    rotation, _ = orthogonal_procrustes(source_centered, target_centered)
    scale = np.linalg.norm(target_centered) / np.linalg.norm(source_centered)
    aligned = source_centered @ rotation * scale + target_mean
    return aligned, {
        "source_mean": source_mean.tolist(),
        "target_mean": target_mean.tolist(),
        "rotation": rotation.tolist(),
        "scale": float(scale),
        "rmse": float(np.sqrt(np.mean((aligned - target) ** 2))),
    }


def orientation_rotation(
    coords: np.ndarray, labels: Sequence[str], target_map: Mapping[str, Tuple[float, float]]
) -> Tuple[np.ndarray, Dict]:
    labels_arr = np.asarray(labels)
    source_centers = []
    target_centers = []
    used = []
    for label, target in target_map.items():
        mask = labels_arr == label
        if mask.sum() == 0:
            continue
        source_centers.append(np.median(coords[mask], axis=0))
        target_centers.append(target)
        used.append(label)
    source_centers_np = np.asarray(source_centers)
    target_centers_np = np.asarray(target_centers)
    source_mean = source_centers_np.mean(axis=0)
    target_mean = target_centers_np.mean(axis=0)
    rotation, _ = orthogonal_procrustes(
        source_centers_np - source_mean, target_centers_np - target_mean
    )
    oriented = (coords - source_mean) @ rotation
    return oriented.astype(np.float32), {
        "labels": used,
        "source_centroid_mean": source_mean.tolist(),
        "target_centroid_mean": target_mean.tolist(),
        "rotation": rotation.tolist(),
        "distance_preserving": True,
    }


def shared_reference_indices(
    event_metadata: pd.DataFrame, reference_meta: pd.DataFrame
) -> np.ndarray:
    ref_lookup = {
        (row.sample_id, int(row.event_index)): idx
        for idx, row in enumerate(reference_meta.itertuples(index=False))
    }
    return np.asarray(
        [ref_lookup[(row.sample_id, int(row.event_index))] for row in event_metadata.itertuples(index=False)],
        dtype=int,
    )




def add_embedding_columns(meta: pd.DataFrame, coords: np.ndarray) -> pd.DataFrame:
    out = meta.copy()
    out["UMAP1"] = coords[:, 0]
    out["UMAP2"] = coords[:, 1]
    out["cell_type_number"] = out["cell_type"].map({x: i + 1 for i, x in enumerate(CELL_TYPES)})
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20260821)
    args = parser.parse_args()
    source_root = SCF_GLOBAL / "Flow_Data" / "fetal flow data"
    panel = SCF_PROJECT_ROOT / "Main_Figure_06/a_M06a"
    plot_data = SCF_SHARED / "Main_Figure_06/a_M06a"
    processed = panel / "03_Source_Data/Transformed_Events"
    plot_data.mkdir(parents=True, exist_ok=True)
    metadata = read_sample_metadata(source_root / "fetal sample information of spectral flow data.xlsx")
    pbmc_records = process_compartment("PBMC", source_root / "23color 20240526 blood.wsp",
                                      source_root / "fcs blood", metadata, processed / "PBMC")
    organ_records = process_compartment("Organ", source_root / "24color 20240526 organ.wsp",
                                       source_root / "fcs organ", metadata, processed / "Organ")
    pbmc_x, pbmc_meta, pbmc_markers = load_cache_set(pbmc_records)
    organ_x, organ_meta, organ_markers = load_cache_set(organ_records)
    pbmc_raw = fit_umap(pbmc_x, args.seed)
    organ_raw = fit_umap(organ_x, args.seed)
    reference = pd.read_csv(SCF_SHARED / "Main_Figure_06/a_M06a/M06a_S14_PBMC_event_coordinates.csv.gz")
    shared_idx = shared_reference_indices(pbmc_meta, reference)
    pbmc_coords, _ = align_similarity(pbmc_raw, reference[["UMAP1", "UMAP2"]].to_numpy()[shared_idx])
    organ_coords, _ = orientation_rotation(organ_raw, organ_meta["cell_type"], ORGAN_TARGETS)
    for name, meta, values, markers, coords in [
        ("PBMC", pbmc_meta, pbmc_x, pbmc_markers, pbmc_coords),
        ("Organ", organ_meta, organ_x, organ_markers, organ_coords),
    ]:
        add_embedding_columns(meta, coords).to_csv(
            plot_data / f"M06a_S14_{name}_event_coordinates.csv.gz", index=False)
        marker_data = pd.DataFrame(values, columns=markers)
        marker_data.insert(0, "event_index", meta["event_index"].to_numpy())
        marker_data.insert(0, "sample_id", meta["sample_id"].to_numpy())
        marker_data.to_csv(SCF_SHARED / "Main_Figure_06/a_M06a" / f"M06a_S14_{name}_transformed_marker_values.csv.gz", index=False)
    pd.DataFrame(pbmc_records + organ_records).drop(columns="cache_path").to_csv(
        SCF_PROJECT_ROOT / "00_Global_Data/Shared_Inputs/Main_Figure_06/a_M06a/M06a_S14_sample_inventory.csv", index=False)
    print("Saved flow UMAP coordinates and sample metadata.")


if __name__ == "__main__":
    main()
