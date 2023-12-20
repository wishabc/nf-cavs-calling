import argparse
import pandas as pd
import scipy.stats as st
from scipy import interpolate
import numpy as np
from tqdm import tqdm


tqdm.pandas()

starting_columns = ['#chr', 'start', 'end', 'ID', 'ref', 'alt']


result_columns = [*starting_columns, 'AAF', 'RAF',
    'mean_BAD', 'nSNPs', 'max_cover', 'mean_cover',
    'footprints_n', 'hotspots_n',
    'group_id',
    'pval_ref_combined',
    'pval_alt_combined',
    'es_combined',
    'logit_es_combined',
    'min_pval',
    'min_fdr'
    ]



def parse_coverage(cov_string):
    try:
        return int(cov_string) if cov_string != 'auto' else 'auto'
    except ValueError:
        print(f'Incorrect coverage threshold provided. {args.max_coverage_tr} not a positive integer or "auto"')
        raise


def logit_es(es, d=1/128):
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

    
def aggregate_pvalues_df(pval_df, groupby_cols=starting_columns):
    pval_df = pval_df.assign(
        **{
            col: pd.NA for col in 
            ['footprints', 'group_id', 'hotspots'] 
            if col not in pval_df.columns
        }
    )
    agg_dict = {
        'nSNPs': ('coverage', 'count'),
        'max_cover': ('coverage', 'max'),
        'hotspots_n': ('hotspots', calc_sum_if_not_minus),
        'footprints_n': ('footprints', calc_sum_if_not_minus),
        'mean_cover': ('coverage', 'mean'),
        'mean_BAD': ('BAD', 'mean'),
        'group_id': ('group_id', 'first'),
        'AAF': ('AAF', 'first'),
        'RAF': ('RAF', 'first')
    }
    for col in groupby_cols:
        if col in agg_dict:
            del agg_dict[col]

    snp_stats = pval_df.groupby(groupby_cols, group_keys=True).agg(**agg_dict)

    agg_pvals = pval_df[[*groupby_cols, 'BAD', 'es', 
        'pval_ref', 'pval_alt', 'inverse_mse', 'coverage']].groupby(
        groupby_cols, group_keys=True
    ).progress_apply(aggregate_pvals)
    return snp_stats.join(agg_pvals).reset_index()


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

def get_min_pval(df, cover_tr, cover_col, pval_cols):
    min_pval = df[pval_cols].min(axis=1) * 2
    ind = (df[cover_col] < cover_tr)
    min_pval[ind] = np.nan
    return min_pval.to_numpy()

def main(pval_df, chrom=None, max_cover_tr=10):
    if chrom is not None:
        pval_df = pval_df[pval_df['#chr'] == args.chrom]
    if pval_df.empty:
        return pd.DataFrame([], columns=result_columns)
    aggr_df = aggregate_pvalues_df(pval_df)
    aggr_df['min_pval'] = get_min_pval(
        aggr_df, 
        cover_tr=max_cover_tr, 
        cover_col='max_cover',
        pval_cols=["pval_ref_combined", "pval_alt_combined"]
    )
    aggr_df['logit_es_combined'] = logit_es(aggr_df['es_combined'])
    aggr_df['min_fdr'] = calc_fdr_pd(aggr_df['min_pval'])
    return aggr_df[result_columns]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('--chrom', help='Chromosome (for parallel execution)', default=None)
    parser.add_argument('--max_coverage_tr', type=str, help="""Threshold for the maximum
                            of coverages of variants aggregated at the same genomic position.
                            Expected to be "auto" or a positive integer""", default='auto')
    args = parser.parse_args()

    coverage_tr = parse_coverage(args.max_coverage_tr)

    pval_df = pd.read_table(args.I, low_memory=False)
    main(pval_df, args.chrom, max_cover_tr=coverage_tr).to_csv(args.O, sep='\t', index=False)
