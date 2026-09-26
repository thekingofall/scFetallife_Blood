import sys
from pathlib import Path
if len(sys.argv) == 1:
    sys.argv.extend([str(Path(__file__).resolve().parents[3] / "00_Global_Data/Fetal_Immune_Atlas_BCR.h5ad"), str(Path(__file__).resolve().parents[1] / "03_Source_Data/BCR_cell_metadata.csv")])
import h5py
import pandas as pd
import numpy as np

def decode(values):
    return [v.decode() if isinstance(v, bytes) else v for v in values]

with h5py.File(sys.argv[1], 'r') as f:
    obs = f['obs']
    index_name = obs.attrs.get('_index', '_index')
    if isinstance(index_name, bytes):
        index_name = index_name.decode()
    index = decode(obs[index_name][:])
    fields = [k for k in obs if k in ('MainID', 'Gestational_Age_Weeks', 'Name', 'Main_Organ')
              or ('call' in k and ('VDJ' in k or 'IGH' in k))]
    data = {}
    for name in fields:
        node = obs[name]
        if isinstance(node, h5py.Group) and 'categories' in node:
            categories = decode(node['categories'][:])
            data[name] = [categories[c] if c >= 0 else None for c in node['codes'][:]]
        elif isinstance(node, h5py.Dataset):
            data[name] = decode(node[:])
    frame = pd.DataFrame(data, index=index)
    frame['Post_Conception_Age_Weeks'] = frame['MainID'].astype(str).str.extract('^[BLTS]([0-9]+(?:\\.[0-9]+)?)(?:_P|\\.P)[0-9]+$', expand=False).astype(float)
    frame.index.name = 'Cellname'
    Path(sys.argv[2]).parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(sys.argv[2])
    print('rows', len(frame), 'fields', list(frame.columns))
