import sys
import pandas as pd


paths = sys.argv[3:]
group_key = sys.argv[2]

with open(sys.argv[1], 'w') as outfile:
    dfs = []
    for fname in paths:
        dfs.append(pd.read_table(fname))
    df = pd.concat(dfs)
    df['group_id'] = group_key
