import argparse
from aggregation import aggregate_pvalues_df, get_min_pval, starting_columns, calc_fdr_pd
import numpy as np
import pandas as pd


def main(tested, pvals, max_cover_tr=15, differential_fdr_tr=0.05, differential_es_tr=0.15):
    constitutive_df = aggregate_pvalues_df(tested)
    constitutive_df['min_pval'] = get_min_pval(
        constitutive_df, 
        cover_tr=max_cover_tr, 
        cover_col='max_cover',
        pval_cols=["pval_ref_combined", "pval_alt_combined"]
    )
    constitutive_df['min_fdr_overall'] = calc_fdr_pd(constitutive_df['min_pval'])

    pvals['differential_fdr'] = calc_fdr_pd(np.exp(pvals['log_p_differential']))
    pvals['differential_es'] = pvals.eval('(F2 + 1) / (F2 + (N - n_groups + 1) / (n_groups - 1))')

    tested_length = len(tested.index)
    tested = tested.merge(
        pvals
    )
    assert len(tested.index) == tested_length

    # set default inividual fdr and find differential snps

    differential_cavs = tested[
        tested.eval(f'differential_fdr <= {differential_fdr_tr} & differential_es >= {differential_es_tr}')
    ].groupby(
        'group_id',
        as_index=False
    ).apply(
        aggregate_pvalues_df
    )
    differential_cavs['min_pval_group'] = get_min_pval(
        differential_cavs, 
        cover_tr=max_cover_tr, 
        cover_col='max_cover',
        pval_cols=["pval_ref_combined", "pval_alt_combined"]
    )
    differential_cavs['min_fdr_group'] = calc_fdr_pd(differential_cavs['min_pval_group'])

    # Group-wise aggregation
    return pvals.merge(
        constitutive_df[[*starting_columns, 'min_pval', 'min_fdr_overall']]
    ).merge(
        differential_cavs[[*starting_columns, 'group_id', 'min_pval_group', 'min_fdr_group']], 
        how='left'
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate ANOVA for tested CAVs')
    parser.add_argument('tested_variants', help='Tested variants')
    parser.add_argument('pvals', help='File with pvals calculated in LRT script')
    parser.add_argument('outpath', help='Outpath to save output to')
    parser.add_argument('--fdr', type=float, help='FDR threshold for differential CAVs', default=0.05)
    parser.add_argument('--es', type=float, help='FDR threshold for differential CAVs', default=0.15)
    
    args = parser.parse_args()
    tested = pd.read_table(args.tested_variants, low_memory=False)
    pvals = pd.read_table(args.pvals)
    res_df = main(tested, pvals, differential_fdr_tr=args.fdr, differential_es_tr=args.es)
    print(len(res_df.index), len(pvals.index))
    res_df.to_csv(args.outpath, sep='\t', index=False)
