import os
import pandas as pd
import argparse
import numpy as np
from scipy.optimize import minimize
from scipy.stats import binom, chi2
from statsmodels.stats.multitest import multipletests
from aggregation import aggregate_pvalues_df, calc_fdr
from tqdm import tqdm
from helpers import starting_columns


tqdm.pandas()

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
    variant_id = df['variant_id'].iloc[0]
    m = len(df['group_id'].unique())
    x = df['x'].to_numpy()
    n = df['n'].to_numpy()
    L0 = censored_binomial_likelihood(x, n, 0.5).sum()
    e1, L1 = get_ml_es_estimation(x, n)
    groups, ess2, stds, DLs, ps_individual, L2 = test_each_group(df)
    p_overall = chi2.logsf(L1 - L0, 1)
    p_differential = chi2.logsf(L2 - L1, m - 1)
    p_zero_int = chi2.logsf(L2 - L0, m)
    return [variant_id, groups, m, e1, ess2, stds, L0, L1 - L0, L2 - L1, ps_individual, p_overall, p_differential, p_zero_int]
    

def censored_binomial_likelihood(xs, ns, p, tr=5):
    b = binom(ns, p)
    return b.logpmf(xs) - np.log(b.cdf(ns - tr) - b.cdf(tr - 1))

def get_ml_es_estimation(x, n):
    def target(p):
        return -censored_binomial_likelihood(x, n, p).sum()
    p = minimize(target, 0.5, bounds=((0.001, 0.999),)).x[0]
    e = np.log2(p) - np.log2(1 - p)
    return e, -target(p)
    

def read_non_aggregated_files(dir_path):
    tables = []
    pval_files = list(os.listdir(dir_path))
    if len(pval_files) == 0:
        print(f'No p-value files found in {dir_path}!')
        raise ValueError
    for file in tqdm(pval_files):
        group_id = file.split('.')[0]
        tb_df = pd.read_table(os.path.join(dir_path, file))
        tb_df['group_id'] = group_id
        tables.append(tb_df)
    res = pd.concat(tables)
    res['variant_id'] = res['#chr'] + '_' + res['end'].astype(str) + '_' + res['alt']
    return res


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


def main(melt_path, min_samples=3, min_groups=2, cover_tr=20):
    melt = read_non_aggregated_files(melt_path)
    melt = melt[melt['#chr'] == 'chr1']
    melt['n'] = melt.eval('ref_counts + alt_counts')
    melt = melt[melt.eval(f'n >= {cover_tr}')]
    melt['x'] = np.round(
        np.where(
            melt['BAD'] == 1, 
            melt['ref_counts'],
            ((1 - 1 / (1 + np.power(2, melt['es']))) * melt['n'])
        )
    )

    testable_pairs = find_testable_pairs(melt, min_samples=min_samples,
        min_groups_per_variant=min_groups)
    # filter only testable variants + cell_types
    tested_melt = melt.merge(
        testable_pairs, on=['variant_id', 'group_id']
    )
    print(f'Testing {len(tested_melt.variant_id.unique())} variants')
    # Total aggregation (find constitutive CAVs)
    constitutive_df = calc_fdr(
        aggregate_pvalues_df(tested_melt, jobs=1, cover_tr=cover_tr)
    ).rename(
        columns={'min_fdr': 'min_fdr_overall'}
    )
    # merge with tested variants
    tested_melt = tested_melt.merge(constitutive_df[['variant_id', 'min_fdr_overall']], how='left')
    print(len(tested_melt.index))

    # LRT (ANOVA-like)
    gb = tested_melt.groupby('variant_id')
    rows = []
    print(len(list(gb.groups)))
    ### TODO: make in parallel
    for g_id in tqdm(list(gb.groups)):
        rows.append(test_snp(gb.get_group(g_id)))
    result = pd.DataFrame(rows,
                     columns=[
                         'variant_id',
                         'group_id',
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
                     ])
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
            lambda x: aggregate_pvalues_df(x, jobs=1, cover_tr=cover_tr)
        )).rename(columns={'min_fdr': 'min_fdr_group'}).reset_index()[['variant_id', 'group_id', 'min_fdr_group']]
    result = result.merge(group_wise_aggregation, how='left')[
        [*starting_columns,
            'group_id', 
            'min_fdr_group',
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
    ]
    return tested_melt, result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate ANOVA for tested CAVs')
    parser.add_argument('input_files_dir', help='Non-aggregated file with tested CAVs')
    parser.add_argument('prefix', help='Prefix to files to save output files into')
    parser.add_argument('--ct', type=int, help='Cover threshold for fdr', default=20)
    parser.add_argument('--min_samples', type=int, help='Number of samples in each group for the variant', default=3)
    parser.add_argument('--min_groups', type=int, help='Number of groups for the variant', default=2)
    args = parser.parse_args()
    min_samples = args.min_samples # of samples
    min_groups_per_variant = args.min_groups
    fdr_cov_tr = args.ct
    df, result = main(
        args.input_files_dir,
        cover_tr=fdr_cov_tr,
        min_samples=min_samples,
        min_groups=min_groups_per_variant)
    
    df.to_csv(f"{args.prefix}.tested.bed", sep='\t', index=False)
    result.to_csv(f"{args.prefix}.cell_selective.bed", sep='\t', index=False)
