"""Prepare RNA, protein, flow, bulk and receptor views for MOFA."""
from pathlib import Path
import os
import subprocess
import sys
code = Path(__file__).resolve().parents[1] / "00_MOFA_Analysis/01_Code"
rscript = os.environ.get("RSCRIPT", "Rscript")
for name in ["01_prepare_single_cell_views.R", "02_prepare_protein_flow_bulk_views.R", "03_prepare_repertoire_views.py", "04_assemble_mofa_views.R"]:
    executable = sys.executable if name.endswith(".py") else rscript
    subprocess.run([executable, str(code / name)], check=True)
