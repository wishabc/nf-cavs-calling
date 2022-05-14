from .helpers import check_states
import sys

for BAD in check_states(sys.argv[2]):
    concat_pval_dfs[concat_pval_dfs['BAD'] == BAD].groupby(['ref_counts', 'alt_counts']).size().reset_index(name='counts').to_csv(
        os.path.join(fit_dir, 'BAD{:.2f}'.format(BAD), 'stats.tsv'), sep='\t', index=False, header=None
    )

