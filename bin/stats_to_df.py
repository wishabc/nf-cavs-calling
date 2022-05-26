import pandas as pd
from helpers import alleles, check_states
import os
import sys


nb_as_params = []
for file in sys.argv[3:]:
    parsed_name, bad, ext = file.split('.')
    allele = parsed_name[-3:]
    bad = float(bad)
    df = pd.read_table(os.path.join(file))
    df['allele'] = allele
    df['BAD'] = bad
    nb_as_params.append(df)
nb_as_params = pd.concat(nb_as_params).reset_index().rename(columns={'index': 'fix_c'})[
    ['allele', 'fix_c', 'BAD', 'r', 'w']
].to_csv(sys.argv[1], sep='\t', index=False)