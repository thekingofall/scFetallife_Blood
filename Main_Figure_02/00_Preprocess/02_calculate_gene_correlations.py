"""Calculate Hotspot autocorrelations, local gene correlations and gene modules."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "00_Settings"))

import argparse
import hotspot
import scanpy as sc
from preprocess_inputs import GLOBAL
parser = argparse.ArgumentParser()
parser.add_argument("--jobs", type=int, default=8)
args = parser.parse_args()
here = Path(__file__).resolve().parent
output = GLOBAL / "Shared_Inputs/Main_Figure_02/c_M02c"
output.mkdir(parents=True, exist_ok=True)
for stage in ["Early", "Late"]:
    data = sc.read_h5ad(here / "02_Processed" / f"M02_PBMC_{stage}_Hotspot.h5ad")
    hs = hotspot.Hotspot(data, layer_key="counts_csc", model="danb", latent_obsm_key="X_pca", umi_counts_obs_key="total_counts")
    hs.create_knn_graph(weighted_graph=False, n_neighbors=25)
    results = hs.compute_autocorrelations(jobs=args.jobs)
    genes = results.loc[results["FDR"] < 0.00001].sort_values("Z", ascending=False).index
    hs.compute_local_correlations(genes, jobs=args.jobs)
    modules = hs.create_modules(min_gene_threshold=20, core_only=False, fdr_threshold=0.01)
    hs.local_correlation_z.to_csv(output / f"{stage.lower()}_local_correlation_z.csv.gz")
    assigned = modules[modules >= 0]
    assigned.rename("module").rename_axis("gene").to_csv(output / f"{stage.lower()}_modules.csv")
    detected = GLOBAL / "Shared_Inputs/Main_Figure_02/d_M02d"
    detected.mkdir(parents=True, exist_ok=True)
    assigned.rename("Module").rename_axis("Gene").to_csv(detected / f"M02d_{stage}_Detected_Modules.csv")
