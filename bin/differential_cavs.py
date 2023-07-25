import argparse
from aggregation import aggregate_pvalues_df, calc_fdr, starting_columns
from statsmodels.stats.multitest import multipletests
import numpy as np
import pandas as pd

def main(tested, pvals, fdr_tr=0.05):
    constitutive_df = calc_fdr(
        aggregate_pvalues_df(tested)
    ).rename(
        columns={'min_fdr': 'min_fdr_overall'}
    )[[*starting_columns, 'min_fdr_overall']]

    pvals['differential_fdr'] = multipletests(
        np.exp(pvals['log_p_differential']),
        method='fdr_bh'
    )[1]
    tested_length = len(tested.index)
    tested = tested.merge(
        pvals[[*starting_columns, 'group_id', 'differential_fdr']]
    )
    assert len(tested.index) == tested_length

    # set default inividual fdr and find differential snps
    differential_cavs = calc_fdr(
        tested[tested['differential_fdr'] <= fdr_tr].groupby(
            'group_id',
            as_index=False
        ).apply(
            aggregate_pvalues_df
        )
    ).rename(
        columns={'min_fdr': 'min_fdr_group'}
    )[[*starting_columns, 'group_id', 'min_fdr_group']]

    # Group-wise aggregation
    return pvals.merge(constitutive_df).merge(differential_cavs, how='left')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate ANOVA for tested CAVs')
    parser.add_argument('tested_variants', help='Tested variants')
    parser.add_argument('pvals', help='File with pvals calculated in LRT script')
    parser.add_argument('outpath', help='Outpath to save output to')
    parser.add_argument('--fdr_tr', type=float, help='FDR threshold for differential CAVs', default=0.05)
    
    args = parser.parse_args()
    tested = pd.read_table(args.tested_variants, low_memory=False)
    pvals = pd.read_table(args.pvals)
    res_df = main(tested, pvals, args.fdr_tr)
    print(len(res_df.index), len(pvals.index))
    res_df.to_csv(args.outpath, sep='\t', index=False)
