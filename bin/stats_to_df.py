import pandas as pd
from helpers import alleles, check_states
import os
import sys


nb_as_params = []
for file in sys.argv[3:]:
    print(file)
for allele in alleles:
    for bad in check_states(sys.argv[1]):
        df = pd.read_table(os.path.join(sys.argv[1],
         'BAD{:.2f}'.format(bad),
          'NBweights_{}.tsv'.format(allele)))
        df['allele'] = allele
        df['BAD'] = bad
        nb_as_params.append(df)
nb_as_params = pd.concat(nb_as_params).reset_index().rename(columns={'index': 'fix_c'})[
    ['allele', 'fix_c', 'BAD', 'r', 'w']
].to_csv(sys.argv[2], sep='\t', index=False)