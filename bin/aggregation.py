import argparse
import pandas as pd
import scipy.stats as st
from scipy import interpolate
import numpy as np
from scipy.special import logit, expit
from tqdm import tqdm


tqdm.pandas()

alleles = {'ref': 'alt', 'alt': 'ref'}
starting_columns = ['#chr', 'start', 'end', 'ID', 'ref', 'alt']


result_columns = [*starting_columns, 'AAF', 'RAF',
    'mean_BAD', 'nSNPs', 'max_cover', 'mean_cover',
    'footprints_n', 'hotspots_n',
    'group_id',
    'pval_ref_weighted',
    'pval_alt_weighted',
    'es_weighted',
    'min_pval',
    'min_fdr'
    ]


def calc_sum_if_not_minus(df_column):
    non_null_vals = [int(x) for x in df_column.tolist() if not pd.isna(x) and x != '-']
    return sum(non_null_vals) if len(non_null_vals) > 0 else '-' 

def aggregate_pvals(df):
    weights = df['inverse_mse']
    pval_ref_weighted = st.combine_pvalues(df['pval_ref'], method='stouffer', weights=weights)[1]
    pval_alt_weighted = st.combine_pvalues(df['pval_alt'], method='stouffer', weights=weights)[1]
    es_weighted = np.average(expit(df['es'] * np.log(2)), weights=weights)
    return pd.Series(
        [pval_ref_weighted, pval_alt_weighted, es_weighted],
        ["pval_ref_weighted", "pval_alt_weighted", "es_weighted"]
        )
    
def aggregate_pvalues_df(pval_df):
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
        pval_df[[*starting_columns, 'BAD', 'es', 'pval_ref', 'pval_alt', 'inverse_mse', 'coverage']].groupby(
            starting_columns
        ).progress_apply(aggregate_pvals)
    ).reset_index()
    # return snp_stats.reset_index()
    # return pval_df[[*starting_columns, 'es', 'min_pval', 'coverage']].groupby(
    #     starting_columns
    # ).progress_apply(
    #     aggregate_es
    # ).join(snp_stats).reset_index()

def qvalue(pvals, bootstrap=False):
    m, pvals = len(pvals), np.asarray(pvals)
    ind = np.argsort(pvals)
    rev_ind = np.argsort(ind)
    pvals = pvals[ind]
    # Estimate proportion of features that are truly null.
    kappa = np.arange(0.05, 0.96, 0.01)
    pik = np.array([sum(pvals > k) / (m*(1-k)) for k in kappa])

    if bootstrap:
        minpi0 = np.quantile(pik, 0.1)
        W = np.array([(pvals>=l).sum() for l in kappa])
        mse = (W / (np.square(m^2) * np.square(1-kappa))) * (1-(W/m)) + np.square((pik-minpi0))
        pi0 = pik[mse==min(mse)][0]
    else:
        cs = interpolate.UnivariateSpline(kappa, pik, k=3, s=None, ext=0)
        pi0 = float(cs(1.))

    pi0 = min(pi0, 1)
    # Compute the q-values.
    qvals = np.zeros(len(pvals))
    qvals[-1] = pi0*pvals[-1]
    for i in np.arange(m-2, -1, -1):
        qvals[i] = min(pi0*m*pvals[i]/float(i+1), qvals[i+1])
    qvals = qvals[rev_ind]
    return qvals


def calc_fdr_pd(pd_series):
    result = np.full(pd_series.shape[0], np.nan)
    ind = pd_series.notna()
    
    if pd_series[ind].shape[0] > 0:  # check if any non-NA p-values exist
        result[ind] = qvalue(pd_series[ind].to_numpy(), bootstrap=True)
    return result


def main(pval_df, chrom=None, max_cover_tr=15):
    if pval_df.empty:
        return pd.DataFrame([], columns=result_columns)
    pval_df = pval_df[pval_df['is_tested']]
    if chrom is not None:
        pval_df = pval_df[pval_df['#chr'] == args.chrom]
    if pval_df.empty:
        return pd.DataFrame([], columns=result_columns)
    aggr_df = aggregate_pvalues_df(pval_df)
    aggr_df['min_pval'] = aggr_df[["pval_ref_weighted", "pval_alt_weighted"]].min(axis=1) * 2
    ind = aggr_df.eval(f'max_cover <= {max_cover_tr} | min_pval > 1')
    aggr_df.loc[ind, 'min_pval'] = pd.NA
    aggr_df['min_fdr'] = calc_fdr_pd(aggr_df['min_pval'])
    return aggr_df[result_columns]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('--chrom', help='Chromosome (for parallel execution)', default=None)
    parser.add_argument('--weights', help='Weights file', default=None)
    parser.add_argument('--max_coverage_tr', type=str, help="""Coverage threshold for recurrent variants.
                                    Expected to be "auto" or a positive integer""", default='auto')
    args = parser.parse_args()
    try:
        coverage_tr = int(args.max_coverage_tr) if args.max_coverage_tr != 'auto' else 'auto'
    except ValueError:
        print(f'Incorrect coverage threshold provided. {args.max_coverage_tr} not a positive integer or "auto"')
        raise
    pval_df = pd.read_table(args.I, low_memory=False)
    #pval_df = pval_df[pval_df['BAD'] <= pval_df[['ref_counts', 'alt_counts']].max(axis=1)/pval_df[['ref_counts', 'alt_counts']].min(axis=1) ]
    weights = pd.read_table(args.weights)
    pval_df = pval_df.merge(weights, on=['BAD', 'coverage'])
    main(pval_df, args.chrom, coverage_tr).to_csv(args.O, sep='\t', index=False)
