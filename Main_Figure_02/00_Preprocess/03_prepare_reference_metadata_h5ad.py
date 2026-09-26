"""Store the annotated fetal/reference cell metadata and UMAP in an AnnData file."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "00_Settings"))

import json
import anndata as ad
import pandas as pd
from scipy import sparse
from preprocess_inputs import GLOBAL, save_h5ad
here = Path(__file__).resolve().parent
types = json.loads((here / "01_Inputs/Reference_Metadata_Types.json").read_text())
metadata = pd.read_csv(
    here / "01_Inputs/Reference_Cell_Metadata_and_UMAP.csv.gz",
    index_col=0, dtype=types, keep_default_na=False, na_values=[""],
)
umap = metadata.pop("UMAP1").to_numpy(), metadata.pop("UMAP2").to_numpy()
import numpy as np
result = ad.AnnData(X=sparse.csr_matrix((len(metadata), 0)), obs=metadata)
result.obsm["X_umap"] = np.column_stack(umap)
colors = json.loads((here / "01_Inputs/Reference_Colors.json").read_text())
result.uns["Last_cell_type_num_colors"] = np.asarray(colors["Last_cell_type_num_colors"])
save_h5ad(result, GLOBAL / "Shared_Inputs/Supplementary_Figures/Supplementary_Figure_05/01_Cell_Metadata/S05_reference_metadata_umap.h5ad")
