import argparse
import pandas as pd
from helpers import get_field_by_ftype, alleles, starting_columns
from calc_pval import result_columns as pval_columns
from scipy.stats import combine_pvalues
from statsmodels.stats.multitest import multipletests
import numpy as np
import multiprocessing as mp


result_columns = starting_columns + ['mean_BAD', 
    '# of SNPs', 'max_cover'] + [get_field_by_ftype(allele, ftype) 
    for allele in alleles for ftype in ('es-ag', 'pval-ag')]

def aggregate_snp(snp_df):
    pvals = {}
    effect_sizes = {}
    for allele in alleles:
        pvals[allele] = logit_aggregate_pvalues(snp_df[get_field_by_ftype(allele)])
        effect_sizes[allele] = aggregate_es(snp_df[get_field_by_ftype(allele, 'es')], 
                                            snp_df[get_field_by_ftype(allele)])
    mean_BAD = snp_df['BAD'].mean()
    return mean_BAD, pvals, effect_sizes


def aggregate_es(es_array, p_array):
    res = [(x, y) for x, y in zip(es_array, p_array) if y != 1 and not pd.isna(y)]
    if len(res) > 0:
        es, p = zip(*res)
        weights = [-1 * np.log10(x) for x in p]
        es_mean = np.round(np.average(es, weights=weights), 3)
    else:
        es_mean = np.nan
    return es_mean


def logit_aggregate_pvalues(pval_list):
    pvalues = np.array([pvalue for pvalue in pval_list if 1 > pvalue > 0])
    if len(pvalues) == 0:
        return 1
    elif len(pvalues) == 1:
        return pvalues[0]
    return combine_pvalues(pvalues, method='mudholkar_george')[1]

def df_to_group(df):
    return df.groupby(['#chr', 'start', 'alt'], as_index=False)

def aggregate_apply(df):
    new_df = df.copy()
    mean_BAD, pvals, effect_sizes = aggregate_snp(df)
    new_df['mean_BAD'] = mean_BAD
    for allele in alleles:
        new_df[get_field_by_ftype(allele, 'pval-ag')] = pvals[allele]
        new_df[get_field_by_ftype(allele, 'es-ag')] = effect_sizes[allele]
    new_df['# of SNPs'] = len(df.index)
    new_df['max_cover'] = df.eval('ref_counts + alt_counts').max()
    new_df = new_df.drop_duplicates(subset=['#chr', 'start', 'alt'])
    return new_df[result_columns]


def aggregate_subgroup(subgroup):
    return pd.concat([aggregate_apply(x.copy()) for x in subgroup])

def aggregate_pvalues_df(pval_df_path, jobs):
    pval_df = pd.read_table(pval_df_path, comment='#', names=pval_columns)
    if pval_df.empty:
        for column in result_columns:
            if column not in pval_df.columns:
                pval_df[column] = None
        return pval_df[result_columns]

    groups = df_to_group(pval_df)
    groups_list = list(groups.groups)
    j = min(jobs,
        max(1, mp.cpu_count()))
    n = len(groups_list) // j
    subgroups = [[groups.get_group(x) for x in groups_list[i: i+n]] for i in range(0, len(groups_list), n)]
    ctx = mp.get_context('forkserver')

    with ctx.Pool(j) as pool:
        results = [pool.apply_async(aggregate_subgroup, (g, )) for g in subgroups]
        snps = []
        for r in results:
            snps.append(r.get())
    return pd.concat(snps)
    
def calc_fdr(aggr_df, max_cover_tr):
    if aggr_df.empty:
        for allele in alleles:
            aggr_df[f"fdrp_bh_{allele}"] = None
        return aggr_df
    mc_filter_array = np.array(aggr_df['max_cover'] >= max_cover_tr)
    for allele in alleles:
        if sum(mc_filter_array) != 0:
            try:
                _, pval_arr, _, _ = multipletests(aggr_df[mc_filter_array][get_field_by_ftype(allele, 'pval-ag')],
                                                                                 alpha=0.05, method='fdr_bh')
            except TypeError:
                print(aggr_df, aggr_df[mc_filter_array][get_field_by_ftype(allele, 'pval-ag')])
                raise
        else:
            pval_arr = []
        fdr_arr = np.empty(len(aggr_df.index), dtype=np.float128)
        fdr_arr[:] = np.nan
        fdr_arr[mc_filter_array] = pval_arr
        aggr_df[f"fdrp_bh_{allele}"] = fdr_arr
    return aggr_df

def main(input_path, out_path, max_covev_tr=10, jobs=1):
    aggr_df = aggregate_pvalues_df(input_path, jobs)
    fdr_df = calc_fdr(aggr_df, max_cover_tr=max_covev_tr)
    fdr_df.to_csv(out_path, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('--mc', type=int, help='Max cover threshold for fdr', default=10)
    parser.add_argument('--jobs', type=int, help='Number of jobs', default=1)
    args = parser.parse_args()
    main(args.I, args.O, args.mc, args.jobs)
