import argparse
import pandas as pd
import scipy.stats as st
from scipy import interpolate
import numpy as np
from tqdm import tqdm


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


def logit_es(es, d=1/128):
    """
    Logit transformation of effect size
    
    Args:
        es: effect size
        d: pseudocount
    
    Returns:
        logit-transformed effect size
    
    """
    return np.log2(es + d) - np.log2(1 - es + d)


def calc_sum_if_not_minus(df_column):
    non_null_vals = [int(x) for x in df_column.tolist() if not pd.isna(x) and x != '-']
    return sum(non_null_vals) if len(non_null_vals) > 0 else '-' 


def aggregate_effect_size(es, weights):
    return np.average(es, weights=weights)


def aggregate_pvals(df):
    weights = df['inverse_mse']
    pval_ref_combined = st.combine_pvalues(df['pval_ref'], method='stouffer', weights=weights)[1]
    pval_alt_combined = st.combine_pvalues(df['pval_alt'], method='stouffer', weights=weights)[1]
    es_combined = aggregate_effect_size(df['es'], weights=weights)
    return pd.Series(
        [pval_ref_combined, pval_alt_combined, es_combined],
        ["pval_ref_combined", "pval_alt_combined", "es_combined"]
    )

def get_min_pval(result, columns):
    min_pval = result[columns].min(axis=1) * 2
    return np.where(np.isnan(min_pval), np.nan, np.minimum(min_pval, 1)).to_numpy()


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
        'nSNPs': ('coverage', 'count'),
        'max_cover': ('max_cover', 'first'),
        'hotspots_n': ('hotspots', calc_sum_if_not_minus),
        'footprints_n': ('footprints', calc_sum_if_not_minus),
        'peaks_n': ('peaks', calc_sum_if_not_minus),
        'mean_cover': ('coverage', 'mean'),
        'mean_BAD': ('BAD', 'mean'),
        'group_id': ('group_id', 'first'),
        'AAF': ('AAF', 'first'),
        'RAF': ('RAF', 'first'),
        'initial_coverage': ('initial_coverage', 'mean'),
        'mean_inverse_mse': ('inverse_mse', 'mean'),
    }
    for col in groupby_cols:
        if col in agg_dict:
            del agg_dict[col]

    snp_stats = pval_df.groupby(groupby_cols, group_keys=True).agg(**agg_dict)
    snp_stats['mean_FMR'] = 1 - snp_stats.eval('mean_cover / initial_coverage')
    agg_pvals = pval_df[[*groupby_cols, 'BAD', 'es', 'is_tested',
        'pval_ref', 'pval_alt', 'inverse_mse', 'coverage']].query('is_tested').groupby(
        groupby_cols, group_keys=True
    ).progress_apply(aggregate_pvals)
    result = snp_stats.join(agg_pvals, how='left').reset_index()
    result['min_pval'] = get_min_pval(result, ['pval_ref_combined', 'pval_alt_combined'])
    result['logit_es_combined'] = logit_es(result['es_combined'])
    return result


# implementation of Storey method for FDR estimation
def qvalue(pvals, bootstrap=False):
    m, pvals = len(pvals), np.asarray(pvals)
    ind = np.argsort(pvals)
    rev_ind = np.argsort(ind)
    pvals = pvals[ind]
    # Estimate proportion of features that are truly null.
    kappa = np.arange(0.05, 0.96, 0.01)
    pik = np.array([sum(pvals > k) / (m * (1-k)) for k in kappa])

    if bootstrap:
        minpi0 = np.quantile(pik, 0.1)
        W = np.array([(pvals >= l).sum() for l in kappa])
        mse = (W / (np.square(m^2) * np.square(1 - kappa))) * (1 - (W / m)) + np.square((pik - minpi0))
        
        if np.any(np.isnan(mse)) or np.any(np.isinf(mse)):
            # Case 1: mse contains NaN or Inf
            mask = np.isnan(mse)  # This will return a boolean mask where True indicates NaN positions

        else:
            # Case 2: mse contains only finite values
            mask = mse == mse.min()
        pi0 = pik[mask][0]
    else:
        cs = interpolate.UnivariateSpline(kappa, pik, k=3, s=None, ext=0)
        pi0 = float(cs(1.))
    
    pi0 = min(pi0, 1)
    # Compute the q-values.
    qvals = np.zeros(len(pvals))
    qvals[-1] = pi0 * pvals[-1]
    for i in np.arange(m - 2, -1, -1):
        qvals[i] = min(pi0 * m * pvals[i]/float(i+1), qvals[i+1])
    qvals = qvals[rev_ind]
    return qvals


def calc_fdr_pd(pd_series):
    ind = pd_series.notna()
    result = np.full_like(pd_series, np.nan, dtype=np.float64)
    if not ind.any():
        return pd_series # all values are NaN
    result[ind] = qvalue(pd_series[ind].to_numpy(), bootstrap=True)
    return result


def main(pval_df, chrom=None):
    if chrom is not None:
        pval_df = pval_df[pval_df['#chr'] == args.chrom]
    pval_df = filter_pval_df(pval_df, args.max_coverage_tr)
    if pval_df.empty or pval_df['is_tested'].sum() == 0:
        return pd.DataFrame([], columns=result_columns)
    
    aggr_df = aggregate_pvalues_df(pval_df)
    aggr_df['min_fdr'] = calc_fdr_pd(aggr_df['min_pval'])
    return aggr_df[result_columns]


def filter_pval_df(df, max_cover_tr=15):
    df['max_cover'] = df.groupby(starting_columns)['coverage'].transform('max')
    df['is_tested'] = (df['max_cover'] >= max_cover_tr) & (df['hotspots'].astype(str) == '1')
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('--chrom', help='Chromosome (for parallel execution)', default=None)
    parser.add_argument('--max_coverage_tr', type=int, help="""Threshold for the maximum
                            of coverages of variants aggregated at the same genomic position.
                            Expected to be a positive integer""", default=15)
    args = parser.parse_args()

    pval_df = pd.read_table(args.I, low_memory=False)
    main(pval_df, args.chrom).to_csv(args.O, sep='\t', index=False)
