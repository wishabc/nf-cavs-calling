import pandas as pd
import argparse
import numpy as np
from scipy.optimize import minimize
from scipy.stats import binom, chi2
from statsmodels.stats.multitest import multipletests
from aggregation import aggregate_pvalues_df, calc_fdr, starting_columns


def test_each_group(df):
    g = df.groupby('group_id')
    ess = []
    stds = []
    liks = []
    deltals = []
    pvals = []
    grps = []
    alpha = np.log(2)/2
    for g_id in g.groups:
        t = g.get_group(g_id)
        es, lhd =  get_ml_es_estimation(t['x'], t['n'])
        lhd0 = censored_binomial_likelihood(t['x'], t['n'], 0.5).sum()
        std = np.cosh(alpha * es) / (alpha * np.sqrt(t['n'].sum()))
        ess.append(es)
        liks.append(lhd)
        stds.append(std)
        deltals.append(lhd - lhd0)
        pvals.append(chi2.logsf(lhd - lhd0, 1))
        grps.append(t['group_id'].iloc[0])
    return grps, ess, stds, deltals, pvals, np.sum(liks)

def test_snp(df):
    variant_id = df.name
    m = df['group_id'].nunique()
    x = df['x'].to_numpy()
    n = df['n'].to_numpy()
    L0 = censored_binomial_likelihood(x, n, 0.5).sum()
    e1, L1 = get_ml_es_estimation(x, n)
    ## rewrite, too slow
    groups, ess2, stds, DLs, ps_individual, L2 = test_each_group(df)
    p_overall = chi2.logsf(L1 - L0, 1)
    p_differential = chi2.logsf(L2 - L1, m - 1)
    p_zero_int = chi2.logsf(L2 - L0, m)
    return pd.DataFrame(
        data=[
            variant_id, groups, m, e1, 
            ess2, stds, 
            L0, L1 - L0, L2 - L1, ps_individual,
            p_overall, p_differential, p_zero_int
        ],
        columns=[
            'variant_id', 'group_id', 'm', 'e1',
            'group_es', 'group_es_std',
            'L0', 'DL1', 'DL2', 'group_pval',
            'p_overall', 'p_differential', 'p_zero_int'
        ]
    )
    

def censored_binomial_likelihood(xs, ns, p, tr=5):
    b = binom(ns, p)
    return b.logpmf(xs) - np.log(b.cdf(ns - tr) - b.cdf(tr - 1))


def get_ml_es_estimation(x, n):
    def target(p):
        return -censored_binomial_likelihood(x, n, p).sum()
    p = minimize(target, 0.5, bounds=((0.001, 0.999),)).x[0]
    e = np.log2(p) - np.log2(1 - p)
    return e, -target(p)


def find_testable_pairs(df, min_samples, min_groups_per_variant):
    # # of samples for particular group with variant_id present
    samples_num = df.value_counts(['group_id', 'variant_id'])
    # for each group, variant is present in >= 3 samples
    samples_num = samples_num[samples_num >= min_samples]
    
    # of variants for particular group
    groups_per_variant = samples_num.reset_index().value_counts('variant_id')

    # variants present in >=2 groups
    tested_ids = groups_per_variant[groups_per_variant >= min_groups_per_variant].index

    testable_variant_group_pairs = samples_num.reset_index().drop(columns=0)

    return testable_variant_group_pairs[
         testable_variant_group_pairs['variant_id'].isin(tested_ids)
    ]


def main(melt, min_samples=3, min_groups=2):
    melt = melt[melt['is_tested']]
    # FIXME
    melt['variant_id'] = melt['#chr'] + '_' + melt['end'].astype(str) + '_' + melt['alt']
    melt['n'] = melt.eval('ref_counts + alt_counts')

    melt['x'] = np.round(
        np.where(
            melt['BAD'] == 1, 
            melt['ref_counts'],
            ((1 - 1 / (1 + np.power(2, melt['es']))) * melt['n'])
        )
    )

    testable_pairs = find_testable_pairs(melt, min_samples=min_samples,
        min_groups_per_variant=min_groups)
    print(testable_pairs)
    # filter only testable variants + cell_types
    tested_melt = melt.merge(
        testable_pairs, on=['variant_id', 'group_id']
    )
    print(f'Testing {tested_melt["variant_id"].nunique()} variants')
    # Total aggregation (find constitutive CAVs)
    print(tested_melt)
    p = aggregate_pvalues_df(tested_melt)
    print(p)
    constitutive_df = calc_fdr(
        aggregate_pvalues_df(tested_melt)
    ).rename(
        columns={'min_fdr': 'min_fdr_overall'}
    )
    # merge with tested variants
    tested_melt = tested_melt.merge(
        constitutive_df[['variant_id', 'min_fdr_overall']], 
        how='left'
    )
    print(len(tested_melt.index))

    # LRT (ANOVA-like)
    rows = tested_melt.groupby('variant_id').progress_apply(test_snp)
    result = pd.DataFrame(rows,
                     )
    result['differential_FDR'] = multipletests(
        np.exp(result['p_differential']),
        method='fdr_bh'
    )[1]
    print(len(result.index))
    
    ## Check and remove 
    # explode result to get by group significance and merge with tested_melt
    result = tested_melt.merge(
        result.explode(['group_id', 'group_es', 'group_es_std', 'group_pval'],
        ignore_index=True), how='left')
    print(len(result.index), result.columns)

    # set default inividual fdr and find differential snps
    differential_idxs = result['differential_FDR'] <= 0.05
    
    # Group-wise aggregation
    group_wise_aggregation = calc_fdr(
        result[differential_idxs].groupby('group_id').progress_apply(
            aggregate_pvalues_df
        )
    ).rename(
        columns={'min_fdr': 'min_fdr_group'}
    ).reset_index()[['variant_id', 'group_id', 'min_fdr_group']]

    result = result.merge(group_wise_aggregation, how='left')[
        [*starting_columns,
            'group_id', 
            'min_fdr_group',
            'min_fdr_overall',
            'differential_FDR',
            'm',
            'e1',
            'group_es',
            'group_es_std',
            'L0',
            'DL1',
            'DL2',
            'group_pval',
            'p_overall',
            'p_differential',
            'p_zero_int'
         ]
    ].drop_duplicates()
    return tested_melt, result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate ANOVA for tested CAVs')
    parser.add_argument('input_data', help='Non-aggregated file with tested CAVs')
    parser.add_argument('prefix', help='Prefix to files to save output files into')
    parser.add_argument('--min_samples', type=int, help='Number of samples in each group for the variant', default=3)
    parser.add_argument('--min_groups', type=int, help='Number of groups for the variant', default=2)
    args = parser.parse_args()

    input_df = pd.read_table(args.input_data)
    df, result = main(
        input_df,
        min_samples=args.min_samples,
        min_groups=args.min_groups
    )
    
    df.to_csv(f"{args.prefix}.differential_tested.bed", sep='\t', index=False)
    result.to_csv(f"{args.prefix}.differential_pvals.bed", sep='\t', index=False)
