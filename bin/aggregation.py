import argparse
import pandas as pd
from helpers import get_pval_field
from helpers import alleles
from scipy.stats import combine_pvalues
from statsmodels.stats.multitest import multipletests
import numpy as np
import multiprocessing as mp


def aggregate_snp(snp_df):
    pvals = {}
    for allele in alleles:
        pvals[allele] = logit_aggregate_pvalues(snp_df[get_pval_field(allele)])
    mean_BAD = snp_df['BAD'].mean()
    return [mean_BAD] + [pvals[field] for field in alleles]


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
    mean_BAD, pval_ref, pval_alt = aggregate_snp(df)
    new_df['mean_BAD'] = mean_BAD
    new_df['logit_pval_ref'] = pval_ref
    new_df['logit_pval_alt'] = pval_alt
    new_df['# of SNPs'] = len(df.index)
    new_df['max_cover'] = df.eval('ref_counts + alt_counts').max()
    new_df = new_df.drop_duplicates(subset=['#chr', 'start', 'alt'])
    new_df = new_df.drop(columns=['ref_counts', 'alt_counts', 'BAD', 'pval_ref', 'pval_alt'])
    return new_df


def aggregate_subgroup(subgroup):
    return pd.concat([aggregate_apply(x.copy()) for x in subgroup])

def aggregate_pvalues_df(pval_df_path, jobs):
    pval_df = pd.read_table(pval_df_path)
    if pval_df.empty:
        return pval_df
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
        return aggr_df
    mc_filter_array = np.array(aggr_df['max_cover'] >= max_cover_tr)
    for allele in alleles:
        if sum(mc_filter_array) != 0:
            _, pval_arr, _, _ = multipletests(aggr_df[mc_filter_array][f"logit_pval_{allele}"],
                                                                                 alpha=0.05, method='fdr_bh')
        else:
            pval_arr = []
        fdr_arr = np.empty(len(aggr_df.index), dtype=np.float128)
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
