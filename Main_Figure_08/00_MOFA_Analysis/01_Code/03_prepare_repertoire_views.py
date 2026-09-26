"""Combine receptor-gene usage with sample-level CDR3 length frequencies."""

from pathlib import Path

import anndata as ad
import pandas as pd


def main():
    root = Path(__file__).resolve().parents[3]
    global_data = root / "00_Global_Data"
    sources = global_data / "Shared_Inputs" / "Main_Figure_08/00_MOFA/01_Source_Data"
    output = Path(__file__).resolve().parent.parent / "02_Data"
    samples = pd.read_csv(output / "MOFA_Sample_Metadata.csv")["MainID"].tolist()
    genes = pd.read_csv(sources / "Receptor_Gene_Usage.csv")
    genes = genes[genes["MainID"].isin(samples)].rename(
        columns={"MainID": "sample_id", "feature": "variable", "view": "type"}
    )
    genes = genes[(genes["type"] != "BCR") | genes["variable"].str.startswith("IGH")].copy()
    schema = pd.read_csv(sources / "Receptor_CDR3_Features.csv")
    results = [genes[["sample_id", "variable", "value", "type"]]]
    adata = ad.read_h5ad(global_data / "Fetal_Immune_Atlas_TCR.h5ad", backed="r")
    obs = adata.obs[["MainID", "Last_cell_type", "IR_VDJ_1_junction_aa"]].copy()
    adata.file.close()
    obs = obs[obs["MainID"].astype(str).isin(samples)]
    for view, cell_label in (("CD4TCR", "CD4"), ("CD8TCR", "CD8")):
        cells = obs[
            obs["Last_cell_type"].astype(str).str.contains(cell_label, regex=False)
        ]
        lengths = cells["IR_VDJ_1_junction_aa"].astype(str).str.len()
        features = schema[schema["view"] == view]
        counts = pd.crosstab(cells["MainID"].astype(str), lengths).reindex(
            index=samples, columns=features["length"].tolist(), fill_value=0
        )
        total = counts.to_numpy().sum()
        if total == 0:
            raise ValueError(f"No CDR3 sequences found for {view}")
        frequencies = counts / total
        frequencies.columns = features["feature"].tolist()
        frequencies.index.name = "sample_id"
        table = frequencies.reset_index().melt(
            id_vars="sample_id", var_name="variable", value_name="value"
        )
        table["type"] = view
        results.append(table)
    pd.concat(results, ignore_index=True).to_csv(output / "Repertoire_Features.csv.gz", index=False)


if __name__ == "__main__":
    main()
