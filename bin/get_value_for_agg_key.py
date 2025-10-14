import pandas as pd
import sys


meta = pd.read_table(sys.argv[1])

aggregation_key = sys.argv[2]

if aggregation_key in ("all", ""):
    meta['all'] = 'all'
    meta['encoded_value'] = 'all'

elif aggregation_key in meta.columns:
    no_dups = meta[[aggregation_key]].drop_duplicates()
    no_dups['path_safe_value'] = no_dups[aggregation_key].astype('category').cat.codes.astype(str)
    mapping = no_dups.set_index(aggregation_key)['path_safe_value'].to_dict()
    meta['encoded_value'] = meta[aggregation_key].map(mapping)
else:
    raise ValueError(f"Aggregation key {aggregation_key} not found in metadata columns")

meta.rename(columns={aggregation_key: 'value'})[
    ['sample_id', 'encoded_value', 'value']
].to_csv(sys.argv[3], index=False, sep='\t')