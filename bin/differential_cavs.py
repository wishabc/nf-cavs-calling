import argparse
from aggregation import aggregate_pvalues_df, get_min_pval, starting_columns, calc_fdr_pd, logit_es
import numpy as np
import pandas as pd
import scipy.stats as st

def main(tested, pvals, differential_fdr_tr=0.05, aggregation_fdr=0.1):

    pvals['differential_fdr'] = calc_fdr_pd(pvals['p_differential'])
    pvals['cell_selective'] = pvals.eval(f'differential_fdr <= {differential_fdr_tr}')
    
    tested_length = len(tested.index)
    tested = tested.merge(pvals)
    assert len(tested.index) == tested_length, f"Length of tested dataframe changed from {tested_length} to {len(tested.index)}"

    # set default inividual fdr and find differential snps

    pvals['pval'] = np.where(
        pvals['cell_selective'], 
        pvals['Pr(>|t|)'], 
        pd.NA
    )

    pvals['fdr_group'] = calc_fdr_pd(pvals['pval'])

    pvals['significant_group'] = pvals.eval(f'cell_selective & fdr_group <= {aggregation_fdr}')
    pvals['has_significant_group'] = pvals.groupby('variant_id')['significant_group'].transform('any')


    pvals['group_es'] = pvals['group_es'] + 0.5
    pvals['logit_group_es'] = logit_es(pvals['group_es'])

    # Group-wise aggregation
    result = pvals.merge(
        constitutive_df[[*starting_columns, 'min_pval', 'min_fdr_overall']]
    )
    initial_len = len(result.index)
    # result = result.merge(
    #     get_category(result)
    # )
    # assert len(result.index) == initial_len

    return result


def get_category(anova, aggregation_fdr=0.1):
    cpy = anova.copy()
    cpy['logit_group_es'] = np.where(
        cpy['fdr_group'].fillna(2) <= aggregation_fdr,
        cpy['group_es'],
        0
    )

    per_variant = cpy.query(f'cell_selective == True').groupby(starting_columns).agg(
        min_es=('logit_group_es', 'min'),
        max_es=('logit_group_es', 'max'),
    )
    
    per_variant = cpy[[*starting_columns, 'has_significant_group', 'cell_selective', 'min_fdr_overall']].drop_duplicates().set_index(starting_columns).join(per_variant)
    per_variant['concordant'] = per_variant.eval('has_significant_group & (min_es * max_es) >= 0')
    per_variant['overall_imbalanced'] = per_variant.eval(f'min_fdr_overall <= {aggregation_fdr}')

    conditions = [
        ~per_variant['overall_imbalanced'] & ~per_variant['cell_selective'], 
        ~per_variant['cell_selective'],
        ~per_variant['has_significant_group'].fillna(False),
        per_variant['concordant'].fillna(False), # concordant                          
    ]

    choices = [
        'not_imbalanced',
        'not_cell_selective',
        'weak_cell_selective',
        'concordant'
    ]
    
    per_variant['category'] = np.select(conditions, choices, default='discordant')
    return per_variant['category'].reset_index()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate ANOVA for tested CAVs')
    parser.add_argument('tested_variants', help='Tested variants')
    parser.add_argument('pvals', help='File with pvals calculated in LRT script')
    parser.add_argument('outpath', help='Outpath to save output to')
    parser.add_argument('--fdr', type=float, help='FDR threshold for differential CAVs', default=0.05)

    args = parser.parse_args()
    tested = pd.read_table(args.tested_variants, low_memory=False)
    pvals = pd.read_table(args.pvals)

    dropped_na_variants = pvals[pvals['p_differential'].isna()]['variant_id'].unique()
    print(dropped_na_variants)
    tested = tested[~tested['variant_id'].isin(dropped_na_variants)]
    pvals = pvals[~pvals['variant_id'].isin(dropped_na_variants)]

    res_df = main(
        tested, pvals, 
        differential_fdr_tr=args.fdr,
    )
    tested['es'] = tested['es'] + 0.5

    tested.to_csv(f"{args.outpath}.tested.bed", sep='\t', index=False)
    tested[tested['variant_id'].isin(dropped_na_variants)].to_csv(f"{args.outpath}.fit_fail.bed", sep='\t', index=False)
    print(len(res_df.index), len(pvals.index))
    res_df.to_csv(f"{args.outpath}.pvals.bed", sep='\t', index=False)
