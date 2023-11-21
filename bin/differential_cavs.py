import argparse
from aggregation import aggregate_pvalues_df, get_min_pval, starting_columns, calc_fdr_pd, logit_es
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
    ]
    print(differential_cavs.columns)
    differential_cavs = differential_cavs.groupby(
        'group_id',
    ).apply(aggregate_pvalues_df)

    print(differential_cavs.columns)
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


def get_category(anova_results):
    cpy = anova_results.copy()
    cpy['abs_max_group_es'] = cpy.groupby('variant_id')['group_es'].transform(lambda x: np.max(np.abs(logit_es(x))))
    cpy['cell_selective'] = cpy.eval('differential_fdr <= 0.05 & differential_es >= 0.2 & abs_max_group_es >= 0.5')

    result = cpy[cpy['cell_selective']][['variant_id', 'min_fdr_group', 'group_es']].groupby('variant_id').progress_apply(get_concordance)
    result = cpy.loc[:, ['variant_id', 'min_fdr_overall', 'differential_fdr', 'differential_es', 'cell_selective']].drop_duplicates().merge(
        result.reset_index(), how='left')
    
    result['overall_imbalanced'] = result['min_fdr_overall'] <= 0.05
    
    result['category'] = 'discordant'
    result.loc[~result['overall_imbalanced'] & ~result['cell_selective'], 'category'] = 'not_imbalanced'
    result.loc[~result['cell_selective'] & pd.isna(result['category']), 'category'] = 'not_cell_selective'
    result.loc[result['concordant'] & pd.isna(result['category']), 'category'] = 'concordant'
    
    result = cpy.loc[:, ['variant_id']].merge(result.reset_index(), on='variant_id')
    return result['category']


def get_concordance(df, fdr_tr=0.05):
    valid_es = df[df.eval(f'min_fdr_group <= {fdr_tr}').fillna(False)]['group_es']
    is_discordant = (valid_es.shape[0] >= 2) and (valid_es.max() - 0.5) * (valid_es.min() - 0.5) < 0
    df['concordant'] = ~is_discordant
    return df[['concordant']].iloc[0]


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
    tested = tested[tested['#chr'] == 'chr12'].copy()
    pvals = pvals[pvals['#chr'] == 'chr12'].copy()
    res_df = main(tested, pvals, differential_fdr_tr=args.fdr, differential_es_tr=args.es)
    print(len(res_df.index), len(pvals.index))
    res_df.to_csv(args.outpath, sep='\t', index=False)
