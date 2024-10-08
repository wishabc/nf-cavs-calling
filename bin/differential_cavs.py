import argparse
from aggregation import aggregate_pvalues_df, check_if_tested, starting_columns, calc_fdr_pd
from AI_statistics.vectorized_estimators import logit_es
import numpy as np
import pandas as pd


def main(tested, pvals, max_cover_tr=15, differential_fdr_tr=0.05, aggregation_fdr=0.1):
    tested_length = len(tested.index)
    tested = tested.merge(pvals)
    assert len(tested.index) == tested_length, f"Length of tested dataframe changed from {tested_length} to {len(tested.index)}"

    tested = check_if_tested(tested, max_cover_tr)
    constitutive_df = aggregate_pvalues_df(tested, starting_columns)
    constitutive_df['min_fdr_overall'] = calc_fdr_pd(constitutive_df['min_pval'])

    constitutive_df['overall_imbalanced'] = constitutive_df.eval(f'min_fdr_overall <= {aggregation_fdr}')
    pvals['differential_fdr'] = calc_fdr_pd(pvals['p_differential'])
    pvals['cell_selective'] = pvals.eval(f'differential_fdr <= {differential_fdr_tr}')
    pvals['n_groups'] = pvals.groupby('variant_id')['group_id'].transform('nunique')
    

    tested['var_es * mse'] = tested.eval('(es - group_es)**2 * inverse_mse')
    tested['es * mse'] = tested.eval('es * inverse_mse')
    mse_estimates = tested.groupby(['variant_id', 'group_id']).agg(
        var_es_sum=('var_es * mse', 'sum'),
        inverse_mse_sum=('inverse_mse', 'sum'),
        mean_group_mse=('inverse_mse', lambda x: np.mean(1/x))
    ).eval(
        '''
        group_es_var = var_es_sum / inverse_mse_sum
        '''
    ).reset_index()[['variant_id', 'group_id', 'group_es_var', 'mean_group_mse']]

    pvals = pvals.drop(
        columns=['group_es_var', 'mean_group_mse'],
        errors='ignore'
    ).merge(
        mse_estimates
    ).merge(
        constitutive_df[
            [
                *starting_columns, 
                'min_pval', 'es_combined',
                'min_fdr_overall', 'overall_imbalanced'
            ]
        ]
    )
    pvals['logit_es_combined'] = logit_es(pvals['es_combined'])

    # set default inividual fdr and find differential snps

    pvals['group_pval'] = np.where(
        pvals['cell_selective'], 
        pvals['Pr(>|t|)'], 
        pd.NA
    )

    pvals['fdr_group'] = calc_fdr_pd(pvals['group_pval'])

    pvals['signif_es'] = np.where(pvals['fdr_group'] < 0.1, pvals['group_es'], 0.5)
    pvals['signif_max_es'] = pvals.groupby('variant_id')['signif_es'].transform('max')
    pvals['signif_min_es'] = pvals.groupby('variant_id')['signif_es'].transform('min')

    pvals['significant_group'] = pvals.eval(f'cell_selective & fdr_group <= {aggregation_fdr}')
    pvals['has_significant_group'] = pvals.groupby('variant_id')['significant_group'].transform('any')

    pvals['concordant'] = pvals.eval('has_significant_group & ((signif_min_es - 0.5) * (signif_max_es - 0.5)) >= 0')
    pvals['discordant'] = pvals.eval('has_significant_group & ((signif_min_es - 0.5) * (signif_max_es - 0.5)) < 0')

    pvals['logit_group_es'] = logit_es(pvals['group_es'])

    return pvals


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate ANOVA for tested CAVs')
    parser.add_argument('tested_variants', help='Tested variants')
    parser.add_argument('pvals', help='File with pvals calculated in LRT script')
    parser.add_argument('outpath', help='Outpath to save output to')
    parser.add_argument('--coverage_tr', help='Coverage threshold for FDR calculation', type=int, default=15)
    parser.add_argument('--fdr', type=float, help='FDR threshold for differential CAVs', default=0.05)

    args = parser.parse_args()
    tested = pd.read_table(args.tested_variants, low_memory=False)
    pvals = pd.read_table(args.pvals)

    dropped_na_variants = pvals[pvals['p_differential'].isna()]['variant_id'].unique()
    failed_fit_df = tested[tested['variant_id'].isin(dropped_na_variants)]
    print(dropped_na_variants.shape, failed_fit_df.shape)
    tested = tested[~tested['variant_id'].isin(dropped_na_variants)]
    pvals = pvals[~pvals['variant_id'].isin(dropped_na_variants)]
    pvals['group_es'] = pvals['group_es'] + 0.5
    tested['es'] = tested['es'] + 0.5
    res_df = main(
        tested, pvals,
        max_cover_tr=args.coverage_tr,
        differential_fdr_tr=args.fdr,
    )

    tested.to_csv(f"{args.outpath}.tested.bed", sep='\t', index=False)
    failed_fit_df.to_csv(f"{args.outpath}.fit_fail.bed", sep='\t', index=False)
    print(len(res_df.index), len(pvals.index))
    res_df.to_csv(f"{args.outpath}.pvals.bed", sep='\t', index=False)
