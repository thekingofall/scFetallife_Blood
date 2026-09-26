"""Export raw counts, normalized expression and cell metadata for Seurat/CellChat."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "00_Settings"))

import gzip
import scanpy as sc
from scipy.io import mmwrite
from preprocess_inputs import GLOBAL, read_counts
here = Path(__file__).resolve().parent / "02_Processed"
here.mkdir(exist_ok=True)
atlas = sc.read_h5ad(GLOBAL / "Fetal_Immune_Atlas.h5ad")
counts = read_counts(atlas.obs_names, atlas.var_names)
for name, matrix in [("counts", counts.X), ("log_normalized", atlas.X)]:
    with gzip.open(here / f"M03_{name}.mtx.gz", "wb") as out:
        mmwrite(out, matrix.T)
atlas.obs.to_csv(here / "M03_Cell_Metadata.csv", index_label="Cellname")
atlas.var_names.to_series().to_csv(here / "M03_Genes.tsv", sep="\t", index=False, header=False)
