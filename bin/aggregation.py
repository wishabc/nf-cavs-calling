import argparse
import pandas as pd
from helpers import alleles, starting_columns
from scipy.stats import combine_pvalues
from statsmodels.stats.multitest import multipletests
import numpy as np
import multiprocessing as mp
import math


keep_columns = [*starting_columns, 'AAF', 'RAF']
result_columns = keep_columns + ['mean_BAD', 
    'nSNPs', 'max_cover',
    'footprints_n', 'hotspots_n',
    'es_weighted_mean', 'es_mean', 
    'logit_pval_ref', 'logit_pval_alt'
    ]


def calc_sum_if_exists(snp_df, column):
    return snp_df[column].sum() if column in snp_df.columns else 0

def aggregate_snp(snp_df):
    pvals = {}
    effect_sizes = {}
    for allele in alleles:
        pvals[allele] = logit_aggregate_pvalues(snp_df[f'pval_{allele}'])

    effect_sizes = aggregate_es(snp_df['es'],
                                snp_df[['pval_ref', 'pval_alt']].min(axis=1),
                                snp_df['coverage'])
    mean_BAD = snp_df['BAD'].mean()
    return mean_BAD, pvals, effect_sizes, calc_sum_if_exists(snp_df, 'hotspots'), calc_sum_if_exists(snp_df, 'footprints')


def expit(x):                                        
   return 1 / (1 + np.exp(-x))

def logit(x):
    return np.log(x) - np.log(1 - x)

def aggregate_es(es_array, p_array, n_array):
    res = [(x, y, n) for x, y, n in zip(es_array, p_array, n_array)
             if y != 1 and not pd.isna(y) and y != 0]
    if len(res) > 0:
        es, p, n = zip(*res)
        weights = [-1 * np.log10(x) for x in p]
        es_weighted_mean = np.average(es, weights=weights)

        sigmas = expit(np.array(es))
        es_mean = logit(np.average(sigmas, weights=n))
    else:
        es_mean = np.nan
        es_weighted_mean = np.nan
    return es_mean, es_weighted_mean


def logit_aggregate_pvalues(pval_list):
    pvalues = np.array([pvalue for pvalue in pval_list if 1 > pvalue > 0])
    if len(pvalues) == 0:
        return 1
    elif len(pvalues) == 1:
        return pvalues[0]
    return combine_pvalues(pvalues, method='mudholkar_george')[1]

def df_to_group(df):
    return df.groupby(starting_columns, as_index=False)

def aggregate_apply(df):
    new_df = df.loc[:, keep_columns].head(1)
    mean_BAD, pvals, effect_sizes, hotspots_n, footprints_n = aggregate_snp(df)
    es_mean, es_weighted_mean = effect_sizes
    new_df['mean_BAD'] = mean_BAD
    new_df['footprints_n'] = footprints_n
    new_df['hotspots_n'] = hotspots_n
    for allele in alleles:
        new_df[f'logit_pval_{allele}'] = pvals[allele]
    new_df['es_mean'] = es_mean
    new_df['es_weighted_mean'] = es_weighted_mean
    new_df['nSNPs'] = len(df.index)
    new_df['max_cover'] = df.eval('coverage').max()
    for extra_column in ('variant_id', 'group_id'):
        if extra_column in df.columns:
            new_df[extra_column] = df.iloc[0].loc[extra_column]
    return new_df


def aggregate_subgroup(subgroup):
    return pd.concat([aggregate_apply(x) for x in subgroup])

def aggregate_pvalues_df(pval_df, jobs, cover_tr):
    pval_df['coverage'] = pval_df.eval('ref_counts + alt_counts')
    if not pval_df.empty:
        pval_df = pval_df[pval_df.eval(f'coverage >= {cover_tr}')]
    if pval_df.empty:
        for column in result_columns:
            if column not in pval_df.columns:
                pval_df[column] = None
        return pval_df[result_columns]

    groups = df_to_group(pval_df)
    groups_list = list(groups.groups)
    j = min(jobs,
        max(1, mp.cpu_count()))
    n = math.ceil(len(groups_list) / j)
    subgroups = [[groups.get_group(x) for x in groups_list[i: i+n]] for i in range(0, len(groups_list), n)]
    snps = []
    if j > 1:
        ctx = mp.get_context('forkserver')
        with ctx.Pool(j) as pool:
            results = [pool.apply_async(aggregate_subgroup, (g, )) for g in subgroups]
            for r in results:
                snps.append(r.get())
    else:
        for subgroup in subgroups:
            snps.append(aggregate_subgroup(subgroup=subgroup))
    return pd.concat(snps)
    
def calc_fdr(aggr_df):
    for allele in alleles:
        if aggr_df.empty:
            fdr_arr = None
        else:
            _, fdr_arr, _, _ = multipletests(
                aggr_df[f'logit_pval_{allele}'],
                alpha=0.05,
                method='fdr_bh'
            )
        aggr_df[f"fdrp_bh_{allele}"] = fdr_arr
    aggr_df['min_fdr'] = aggr_df[[f'fdrp_bh_{x}' for x in alleles]].min(axis=1)
    return aggr_df

def main(input_path, out_path, cover_tr=10, jobs=1):
    pval_df = pd.read_table(input_path)
    aggr_df = aggregate_pvalues_df(pval_df, jobs, cover_tr)
    fdr_df = calc_fdr(aggr_df)
    fdr_df.to_csv(out_path, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('--ct', type=int, help='Cover threshold for fdr', default=10)
    parser.add_argument('--jobs', type=int, help='Number of jobs', default=1)
    args = parser.parse_args()
    main(args.I, args.O, args.ct, args.jobs)
