"""Score doublets per sample, normalize RNA and assign CellTypist labels."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "00_Settings"))

import anndata as ad
import celltypist
import scanpy as sc
import scrublet as scr
from preprocess_inputs import sample_inputs, read_sample_counts, save_h5ad

output = Path(__file__).resolve().parent / "02_Processed"
blocks = []
for row in sample_inputs().itertuples(index=False):
    data = read_sample_counts(row)
    detector = scr.Scrublet(data.X)
    scores, _ = detector.scrub_doublets(
        min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
    predicted = detector.call_doublets(threshold=0.8)
    data.obs["doublet_scores"] = scores
    data.obs["predicted_doublets"] = predicted
    blocks.append(data)
merged = ad.concat(blocks, join="outer", merge="first", fill_value=0)
save_h5ad(merged, output / "M01_Raw_Counts.h5ad")
sc.pp.normalize_total(merged, target_sum=1e4)
sc.pp.log1p(merged)
model = ROOT / "00_Global_Data/Raw_Inputs/CellTypist/Immune_All_Low.pkl"
prediction = celltypist.annotate(merged, model=str(model), majority_voting=True)
save_h5ad(prediction.to_adata(), output / "M01_CellTypist.h5ad")
