import argparse
import pandas as pd
import scipy.stats as st
from statsmodels.stats.multitest import multipletests
import numpy as np
from scipy.special import logit, expit
from tqdm import tqdm
import json
tqdm.pandas()

alleles = {'ref': 'alt', 'alt': 'ref'}
starting_columns = ['#chr', 'start', 'end', 'ID', 'ref', 'alt']


result_columns = [*starting_columns, 'AAF', 'RAF',
    'mean_BAD', 'nSNPs', 'max_cover', 'mean_cover',
    'footprints_n', 'hotspots_n',
    #'es_weighted_mean', 'es_mean', 
    #'logit_pval_ref', 'logit_pval_alt',
     'group_id',
    #'fdrp_bh_ref', 'fdrp_bh_alt',
    'aggregated_pval_ref',
    'aggregated_pval_alt',
    'aggregated_pval',
     'min_fdr'
    ]


def calc_sum_if_not_minus(df_column):
    non_null_vals = [int(x) for x in df_column.tolist() if not pd.isna(x) and x != '-']
    return sum(non_null_vals) if len(non_null_vals) > 0 else '-' 


def aggregate_es(stat):
    valid_index = ~pd.isna(stat['min_pval']) & (stat['min_pval'] != 1) & (stat['min_pval'] != 0)
    stat = stat[valid_index]
    es_column = stat['es']
    if es_column.empty:
        es_mean = np.nan
        es_weighted_mean = np.nan
    else:
        es_weighted_mean = np.average(es_column.to_numpy(), weights=-np.log10(stat['min_pval']))
        es_mean = logit(np.average(expit(es_column.to_numpy()), weights=stat['coverage']))
    
    return pd.Series([es_mean, es_weighted_mean], ["es_mean", "es_weighted_mean"])


def aggregate_pvalues(pval_list, method='stouffer'):
    pvalues = pval_list
    # pvalues = np.array([pvalue for pvalue in pval_list if 1 > pvalue > 0])
    # if len(pvalues) == 0:
    #     return 1
    # elif len(pvalues) == 1:
    #     return pvalues[0]
    return st.combine_pvalues(pvalues, method=method,)[1]

def aggregate_pvals_stf(df, weights_dict):
    # weights = np.power(df['coverage'], 2)
    weights = df[['BAD', 'coverage']].apply(lambda row: weights_dict[str(float(row['BAD']))][str(int(row['coverage']))], axis=1).to_numpy()
    weights = np.power(weights, 2)
    pval_ref_weighted = st.combine_pvalues(df['pval_ref'], method='stouffer', weights=weights)[1]
    pval_alt_weighted = st.combine_pvalues(df['pval_alt'], method='stouffer', weights=weights)[1]
    #pval_weighted = st.combine_pvalues(df['min_pval'].fillna(1), method='stouffer', weights=weights)[1]
    
    # pval_ref_weighted2 = st.combine_pvalues(df['pval_ref'], method='stouffer', weights=weights)[1]
    # pval_alt_weighted2 = st.combine_pvalues(df['pval_alt'], method='stouffer', weights=weights)[1]
    # pval_weighted2 = st.combine_pvalues(df['min_pval'], method='stouffer', weights=weights)[1]
    return pd.Series(
        [pval_ref_weighted, pval_alt_weighted],
        ["pval_ref_weighted", "pval_alt_weighted"]
        )
    
def aggregate_pvalues_df(pval_df, weights):
    pval_df = pval_df.assign(
        **{
            col: pd.NA for col in 
            ['footprints', 'group_id', 'hotspots'] 
            if col not in pval_df.columns
        }
    )
    snp_stats = pval_df.groupby(starting_columns).agg(
        nSNPs=('coverage', 'count'),
        max_cover=('coverage', 'max'),
        hotspots_n=('hotspots', calc_sum_if_not_minus),
        footprints_n=('footprints', calc_sum_if_not_minus),
        mean_cover=('coverage', 'mean'),
        mean_BAD=('BAD', 'mean'),
        group_id=('group_id', 'first'),
        AAF=('AAF', 'first'),
        RAF=('RAF', 'first')
    )
    return snp_stats.join(
        pval_df[[*starting_columns, 'BAD', 'pval_ref', 'pval_alt', 'min_pval', 'coverage']].groupby(
            starting_columns
        ).progress_apply(
            lambda x: aggregate_pvals_stf(x, weights)
        )
    ).reset_index()
    # return snp_stats.reset_index()
    # return pval_df[[*starting_columns, 'es', 'min_pval', 'coverage']].groupby(
    #     starting_columns
    # ).progress_apply(
    #     aggregate_es
    # ).join(snp_stats).reset_index()


def calc_fdr(aggr_df, prefix='aggregated_'):
    aggr_df["min_fdr"] = multipletests(
            aggr_df[f'{prefix}pval'],
            alpha=0.05,
            method='fdr_bh'
        )[1]
    return aggr_df


def main(pval_df, chrom=None, weights=None):
    if pval_df.empty:
        return pd.DataFrame([], columns=result_columns)
    pval_df = pval_df[pval_df['is_tested']]
    if chrom is not None:
        pval_df = pval_df[pval_df['#chr'] == args.chrom]
    if pval_df.empty:
        return pd.DataFrame([], columns=result_columns)
    aggr_df = aggregate_pvalues_df(pval_df, weights)
    #res_df = calc_fdr(aggr_df)
    #return res_df
    return aggr_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('--chrom', help='Chromosome (for parallel execution)', default=None)
    parser.add_argument('--weights', help='Weights file', default=None)
    args = parser.parse_args()
    pval_df = pd.read_table(args.I, low_memory=False)
    with open(args.weights) as f:
        weights = json.load(f)
    main(pval_df, args.chrom, weights).to_csv(args.O, sep='\t', index=False)
