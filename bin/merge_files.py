import sys
import pandas as pd


paths = sys.argv[3:]
group_key = sys.argv[2]

with open(sys.argv[1], 'w') as outfile:
    df = pd.concat([pd.read_table(fname) for fname in paths])
    df['group_id'] = group_key
