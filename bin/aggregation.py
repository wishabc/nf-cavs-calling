import argparse
import pandas as pd
import numpy as np
from tqdm import tqdm

from AI_statistics.vectorized_estimators import logit_es, aggregate_effect_size, stouffer_combine_log_pvals, log_pval_both, qvalue


tqdm.pandas()


starting_columns = ['#chr', 'start', 'end', 'ID', 'ref', 'alt']
result_columns = [
    *starting_columns, 'AAF', 'RAF', 
    'mean_FMR', 'mean_BAD',
    'nSNPs', 'max_cover', 'mean_cover', 'mean_inverse_mse', 
    'footprints_n', 'hotspots_n', 'peaks_n',
    'group_id',
    'pval_ref_combined',
    'pval_alt_combined',
    'es_combined',
    'logit_es_combined',
    'min_pval',
    'min_fdr'
]


def calc_sum_if_not_minus(df_column):
    non_null_vals = [int(x) for x in df_column.tolist() if not pd.isna(x) and x != '-']
    return sum(non_null_vals) if len(non_null_vals) > 0 else '-' 


def aggregate_pvals(df):
    weights = df['inverse_mse']
    log_pval_ref_combined = stouffer_combine_log_pvals(df['pval_ref'], weights)
    log_pval_alt_combined = stouffer_combine_log_pvals(df['pval_alt'], weights)
    log_pval_both_combined = log_pval_both(log_pval_ref_combined, log_pval_alt_combined)
    es_combined = aggregate_effect_size(df['es'], weights=weights)
    return pd.Series(
        [np.exp(log_pval_ref_combined), np.exp(log_pval_alt_combined), es_combined, np.exp(log_pval_both_combined)],
        ["pval_ref_combined", "pval_alt_combined", "es_combined", 'min_pval']
    )


def aggregate_pvalues_df(pval_df, groupby_cols=None):
    if groupby_cols is None:
        groupby_cols = starting_columns

    pval_df = pval_df.assign(
        **{
            col: pd.NA for col in 
            ['footprints', 'group_id', 'hotspots', 'peaks'] 
            if col not in pval_df.columns
        }
    )
    pval_df['initial_coverage'] = pval_df.eval('coverage / (1 - FMR)')
    agg_dict = {
        'AAF': ('AAF', 'first'),
        'RAF': ('RAF', 'first'),

        'nSNPs': ('coverage', 'count'),
        'max_cover': ('coverage', 'max'),
        'hotspots_n': ('hotspots', calc_sum_if_not_minus),
        'footprints_n': ('footprints', calc_sum_if_not_minus),
        'peaks_n': ('peaks', calc_sum_if_not_minus),
        'mean_cover': ('coverage', 'mean'),
        'mean_BAD': ('BAD', 'mean'),
        'group_id': ('group_id', 'first'),
        'initial_coverage': ('initial_coverage', 'mean'),
        'mean_inverse_mse': ('inverse_mse', 'mean'),
    }
    for col in groupby_cols:
        if col in agg_dict:
            del agg_dict[col]

    snp_stats = pval_df.groupby(groupby_cols, group_keys=True).agg(**agg_dict)
    snp_stats['mean_FMR'] = 1 - snp_stats.eval('mean_cover / initial_coverage')
    agg_pvals = pval_df[[*groupby_cols, 'BAD', 'es', 'is_tested',
        'pval_ref', 'pval_alt', 'inverse_mse', 'coverage']].groupby(
        groupby_cols, group_keys=True
    ).progress_apply(aggregate_pvals)
    result = snp_stats.join(agg_pvals, how='left').reset_index()
    result['logit_es_combined'] = logit_es(result['es_combined'])
    return result


def calc_fdr_pd(pd_series):
    ind = pd_series.notna()
    result = np.full_like(pd_series, np.nan, dtype=np.float64)
    if not ind.any():
        return pd_series # all values are NaN
    result[ind] = qvalue(pd_series[ind].to_numpy(), bootstrap=True)
    return result


def main(pval_df, max_cover_tr=15, chrom=None):
    if chrom is not None:
        pval_df = pval_df[pval_df['#chr'] == args.chrom]
    pval_df = check_if_tested(pval_df, max_cover_tr=max_cover_tr).query('is_tested').reset_index(drop=True)
    if pval_df.empty:
        return pd.DataFrame([], columns=result_columns)
    
    aggr_df = aggregate_pvalues_df(pval_df)
    aggr_df['min_fdr'] = calc_fdr_pd(aggr_df['min_pval'])
    return aggr_df[result_columns]


def check_if_tested(df, max_cover_tr=15):
    if "hotspots" in df.columns:
        df['is_tested'] = df['hotspots'].astype(str).isin(['1', '-'])
    else:
        df['is_tested'] = True
    df['max_cover'] = df.query('is_tested').groupby(starting_columns)['coverage'].transform('max')
    df['is_tested'] = df.eval(f'max_cover >= {max_cover_tr} & is_tested')
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('--chrom', help='Chromosome (for parallel execution)', default=None)
    parser.add_argument('--max_coverage_tr', type=int, help="""Threshold for the highest coverage of the variants aggregated at the same genomic position.
    Expected to be a positive integer""", default=15)
    args = parser.parse_args()

    pval_df = pd.read_table(args.I, low_memory=False)
    main(pval_df, args.max_coverage_tr, args.chrom).to_csv(args.O, sep='\t', index=False)
