"""Define exact and sequence-similarity TCR clonotypes and calculate clone sizes."""
from pathlib import Path
import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "00_Settings"))

import scanpy as sc
import scirpy as ir
from preprocess_inputs import GLOBAL, save_h5ad
here = Path(__file__).resolve().parent
data = sc.read_h5ad(here / "02_Processed/M04_Paired_TCR.h5ad")
ir.pp.ir_dist(data, sequence="aa")
ir.tl.define_clonotypes(data, receptor_arms="all", dual_ir="primary_only")
ir.tl.clonotype_network(data, min_cells=2, random_state=42)
ir.pp.ir_dist(data, metric="alignment", sequence="aa", cutoff=15)
ir.tl.define_clonotype_clusters(data, sequence="aa", metric="alignment", receptor_arms="all", dual_ir="any")
ir.tl.define_clonotype_clusters(data, sequence="aa", metric="alignment", receptor_arms="all", dual_ir="any", same_v_gene=True, key_added="cc_aa_alignment_same_v")
ir.tl.clonal_expansion(data)
save_h5ad(data, GLOBAL / "Fetal_Immune_Atlas_TCR.h5ad")
