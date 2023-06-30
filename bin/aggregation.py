import argparse
import pandas as pd
from helpers import alleles, starting_columns
from scipy.stats import combine_pvalues
from statsmodels.stats.multitest import multipletests
from calc_pval_binom import calc_pval
import numpy as np
from tqdm import tqdm

tqdm.pandas()

keep_columns = [*starting_columns, 'AAF', 'RAF']
result_columns = keep_columns + [
    'mean_BAD', 'nSNPs', 'max_cover', 'mean_cover',
    'footprints_n', 'hotspots_n',
    'es_weighted_mean', 'es_mean', 
    'logit_pval_ref', 'logit_pval_alt', 'group_id'
    ]

def calc_min_cover_by_BAD(BAD, es=1, pvalue_tr=0.05, allele_tr=5, cmax=1000):
    covs = np.arange(allele_tr * 2, cmax + 1)
    x = covs / (1 + np.power(2.0, -es) / BAD)
    mask = (allele_tr <= x) & (x <= covs - allele_tr)
    x = x[mask]
    covs = covs[mask]
    BADs = np.full(x.shape, BAD)
    _, _, sigs, _ = calc_pval(covs, x, BADs, allele_tr=allele_tr, modify_w=True, smooth=True)
    return np.min(covs[sigs < pvalue_tr])


def calc_sum_if_not_minus(df_column):
    non_null_vals = [int(x) for x in df_column.tolist() if not pd.isna(x) and x != '-']
    return sum(non_null_vals) if len(non_null_vals) > 0 else '-' 


def expit(x):                                        
   return 1 / (1 + np.exp(-x))

def logit(x):
    return np.log(x) - np.log(1 - x)

def aggregate_es(df, pval_df):
    print(df.index)
    es = pval_df[df.index]
    df = df[~pd.isna(df['min_pval']) & (df['min_pval'] != 1) & (df['min_pval'] != 0)]

    if df.empty:
        es_mean = np.nan
        es_weighted_mean = np.nan
    else:
        es_weighted_mean = np.average(df['es'], weights=-np.log10(df['min_pval']))
        es_mean = logit(np.average(expit(df['es']), weights=df['coverage']))
    
    return pd.Series([es_mean, es_weighted_mean], ['es_mean', 'es_weighted_mean'])

def logit_aggregate_pvalues(pval_list):
    pvalues = np.array([pvalue for pvalue in pval_list if 1 > pvalue > 0])
    if len(pvalues) == 0:
        return 1
    elif len(pvalues) == 1:
        return pvalues[0]
    return combine_pvalues(pvalues, method='mudholkar_george')[1]

def df_to_group(df):
    return df.groupby(keep_columns)

def flatten_colname(data):
    if isinstance(data, tuple):
            return '_'.join(flatten_colname(val) for val in data)
    else:
        return data

def aggregate_pvalues_df(pval_df):
    groups = df_to_group(pval_df)
    snp_stats = groups.agg(
        max_cover=('coverage', 'max'),
        logit_pval_alt=('pval_alt', logit_aggregate_pvalues),
        logit_pval_ref=('pval_ref', logit_aggregate_pvalues),
        hotspots_n=('hotspots', calc_sum_if_not_minus),
        footprints_n=('footprints', calc_sum_if_not_minus),
        mean_cover=('coverage', 'mean'),
        mean_BAD=('BAD', 'mean'),
        group_id=('group_id', 'first'),
        es=('es', lambda x: aggregate_es(x, pval_df))
    )
    t = groups.apply(
            aggregate_es, pval_df
        ).join(
            snp_stats
        ).reset_index()
    t.columns = [flatten_colname(col) for col in t.columns.values]
    print(t.columns)
    return t.rename(columns={
            'coverage_max': 'max_cover',
            'coverage_mean': 'mean_cover',
            'ref_counts_count': 'nSNPs',
            'BAD_mean': 'mean_BAD',
            'footprints_calc_sum_if_not_minus': 'footrpints_n',
            'hotspots_calc_sum_if_not_minus': 'hotspots_n',
            'pval_ref_logit_aggregate_pvalues': 'logit_pval_ref',
            'pval_alt_logit_aggregate_pvalues': 'logit_pval_alt',
            'group_id_first': 'group_id'
        }
    )
    
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

def main(input_path, coverage_tr):
    pval_df = pd.read_table(input_path)
    if pval_df.empty:
        for column in result_columns:
            if column not in pval_df.columns:
                pval_df[column] = None
        return pval_df[result_columns]
    for column in ('variant_id', 'group_id', 'hotspots', 'footprints'):
        if column not in pval_df.columns:
            pval_df[column] = pd.NA

    pval_df['coverage'] = pval_df.eval('ref_counts + alt_counts')
    pval_df['min_pval'] = pval_df[['pval_ref', 'pval_alt']].min(axis=1)

    if coverage_tr == 'auto':
        by_BAD_coverage_tr = {x: calc_min_cover_by_BAD(x) for x in pval_df['BAD'].unique()}
        pval_df = pval_df[pval_df['coverage'] >= pval_df['BAD'].apply(lambda x: by_BAD_coverage_tr[x])]
    else:    
        pval_df = pval_df[pval_df.eval(f'coverage >= {coverage_tr}')]

    
    aggr_df = aggregate_pvalues_df(pval_df)
    return calc_fdr(aggr_df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('--ct', type=str, help="""Coverage threshold for individual variants to be aggregated.
                                                Should be "auto" or positive integer""", default='auto')
    args = parser.parse_args()
    try:
        coverage_tr = int(args.ct) if args.ct != 'auto' else 'auto'
    except ValueError:
        print(f'Incorrect coverage threshold provided. {args.ct} not a positive integer or "auto"')
        raise
    final_df = main(args.I, coverage_tr)
    final_df.to_csv(args.O, sep='\t', index=False)
