import argparse
from aggregation import aggregate_pvalues_df, get_min_pval, starting_columns, calc_fdr_pd, logit_es
import numpy as np
import pandas as pd
import scipy.stats as st

def main(tested, pvals, max_cover_tr=15, differential_fdr_tr=0.05, differential_es_tr=0.15):
    constitutive_df = aggregate_pvalues_df(tested)
    constitutive_df['min_pval'] = get_min_pval(
        constitutive_df, 
        cover_tr=max_cover_tr, 
        cover_col='max_cover',
        pval_cols=["pval_ref_combined", "pval_alt_combined"]
    )
    constitutive_df['min_fdr_overall'] = calc_fdr_pd(constitutive_df['min_pval'])

    pvals['differential_fdr'] = calc_fdr_pd(pvals['p_differential'])
    pvals['differential_es'] = pvals.eval('var_group_id / (var_residuals + var_group_id + var_indiv_id)')

    tested_length = len(tested.index)
    tested = tested.merge(pvals)
    assert len(tested.index) == tested_length

    # set default inividual fdr and find differential snps

    differential_cavs = tested.query(
        f'differential_fdr <= {differential_fdr_tr} 
        & differential_es >= {differential_es_tr}'
    )

    differential_cavs = aggregate_pvalues_df(
        differential_cavs, 
        groupby_cols=[*starting_columns, 'group_id']
    )

    differential_cavs['pval_group'] = st.t.sf(
        abs(differential_cavs.eval('(group_es - 0.5)/group_es_std**2')),
        df=differential_cavs['samples_count']
    )

    differential_cavs['fdr_group'] = calc_fdr_pd(differential_cavs['pval_group'])

    # Group-wise aggregation
    result = pvals.merge(
        constitutive_df[[*starting_columns, 'min_pval', 'min_fdr_overall']]
    )

    initial_len = len(result.index)
    result = result.merge(
        get_category(result, differential_fdr_tr=0.05)
    )
    assert len(result.index) == initial_len

    return result


def get_category(anova_results, differential_fdr_tr=0.05, aggregation_fdr=0.1):
    cpy = anova_results.copy()

    cpy['logit_group_es'] = logit_es(cpy['group_es'])
    cpy['logit_group_es'] = np.where(
        cpy['fdr_group'].notna(),
        cpy['group_es'],
        0
    )
    cpy['abs_logit_es'] = np.abs(cpy['logit_group_es'])

    cpy['cell_selective'] = cpy.eval(f'differential_fdr <= {differential_fdr_tr}')
    
    cpy['significant_group'] = cpy.eval(f'cell_selective & fdr_group <= {aggregation_fdr}')

    result = cpy.query(f'cell_selective == True').groupby(starting_columns).agg(
        strong_cell_selective=('significant_group', 'any'),
        cell_selective=('cell_selective', 'any'),
        min_es=('logit_group_es', 'min'),
        max_es=('logit_group_es', 'max'),
    )
    result['concordant'] = result.eval('strong_cell_selective & min_es * max_es >= 0')
    result = cpy[[*starting_columns, 'cell_selective', 'min_fdr_overall', 'overall_es']].drop_duplicates(
        ).set_index(
            starting_columns
        ).join(result)
    result['overall_imbalanced'] = result.eval(f'min_fdr_overall <= {aggregation_fdr}')
    
    conditions = [
        ~result['overall_imbalanced'] & ~result['cell_selective'], # not_imbalanced
        ~result['cell_selective'],                                 # not_cell_selective
        ~result['strong_cell_selective'],
        result['concordant'].fillna(True)                          # concordant
    ]

    choices = [
        'not_imbalanced',
        'not_cell_selective',
        'weak_cell_selective',
        'concordant'
    ]
    
    result['category'] = np.select(conditions, choices, default='discordant')
    return result['category'].reset_index()


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

    dropped_na_variants = pvals[pvals['p_differential'].isna()]['variant_id'].unique()
    
    tested = tested[~tested['variant_id'].isin(dropped_na_variants)]
    pvals = pvals[pvals['p_differential'].notna()]


    res_df = main(
        tested, pvals, 
        differential_fdr_tr=args.fdr,
        differential_es_tr=args.es
    )
    print(len(res_df.index), len(pvals.index))
    res_df.to_csv(args.outpath, sep='\t', index=False)
