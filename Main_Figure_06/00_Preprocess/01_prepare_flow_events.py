"""Apply the stored FlowJo gates, transform marker values and calculate flow UMAP coordinates."""
from pathlib import Path
import runpy
script = Path(__file__).resolve().parents[1] / "a_M06a/01_Code/01_gate_flow_events_and_compute_umap.py"
runpy.run_path(str(script), run_name="__main__")
