from helpers import check_states
import argparse
import pandas as pd
import os

def main(badmaps, fit_dir):
    concat_pval_dfs = pd.read_table(badmaps)
    if concat_pval_dfs.empty:
        raise ValueError(f'No data in {badmaps}')
    for BAD in concat_pval_dfs['BAD'].unique():
        bad_dir = os.path.join(fit_dir, 'BAD{:.2f}'.format(BAD))
        if not os.path.exists(bad_dir):
            os.mkdir(bad_dir)
        concat_pval_dfs[concat_pval_dfs['BAD'] == BAD].groupby(['ref_counts', 'alt_counts']).size().reset_index(name='counts').to_csv(
            os.path.join(bad_dir, 'stats.tsv'), sep='\t', index=False, header=None
        )
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collect NB stats')
    parser.add_argument('-O', help='Dir to save calculated p-value into')
    parser.add_argument('-b', help='Badmaps files')
    args = parser.parse_args()
    main(args.b, args.O)
