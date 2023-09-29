import pandas as pd
from scipy.stats import binom
from scipy.special import betainc
import argparse
import numpy as np
from scipy.special import expit
from statsmodels.stats.multitest import multipletests

from tqdm import tqdm

tqdm.pandas()


updated_columns = ['w', 'es', 'pval_ref', 'pval_alt', 'is_tested']

class CalcImbalance:
    def __init__(self, allele_tr, modify_w):
        self.allele_tr = allele_tr
        self.modify_w = modify_w
    
    def calc_pval(self, n, k, BADs, smooth=False):
        p = 1 / (BADs + 1)
        ws = np.full(k.shape, 0.5, dtype=np.float64)
        if self.modify_w:
            ws = self.modify_w_binom(k, n, p, ws)
        
        es = self.odds_es(k, n, BADs, ws)
        pval_ref = self.censored_binom_pvalue(k, n, p, ws, smooth=smooth)
        pval_alt = self.censored_binom_pvalue(n - k, n, p, 1 - ws, smooth=smooth)
        return [ws, es, pval_ref, pval_alt]

    def calc_min_cover_by_BAD(self, BAD, pvalue_tr=0.05, cmax=1000):
        covs = np.arange(self.allele_tr * 2, cmax + 1)
        x = covs - self.allele_tr
        BADs = np.full(covs.shape, BAD)
        sigs = self.calc_pval(covs, x, BADs, smooth=True)[2]
        return np.min(covs[sigs <= pvalue_tr])

    @staticmethod
    def cdf(x, n, p, smooth):
        if smooth:
            return betainc(n - x, x + 1, 1 - p)
        else:
            return binom.cdf(x, n, p)

    @staticmethod
    def sf(x, n, p, smooth):
        if smooth:
            return betainc(x + 1, n - x, p)        
        else:
            return binom.sf(x, n, p)


    def binom_sf_mixture(self, x, n, p, weight, smooth):
        return (1 - weight) * self.sf(x, n, p, smooth) + weight * self.sf(x, n, 1 - p, smooth)

    def binom_cdf_mixture(self, x, n, p, w, smooth):
        return (1 - w) * self.cdf(x, n, p, smooth) + w * self.cdf(x, n, 1 - p, smooth)


    def censored_binom_pvalue(self, x, n, p, w, smooth):
        if self.allele_tr > 0:
            norm_coef = self.binom_cdf_mixture(n - (self.allele_tr - 1), n, p, w, smooth) - self.binom_cdf_mixture(self.allele_tr - 1, n, p, w, smooth)
            return (self.binom_sf_mixture(x - 1, n, p, w, smooth) - self.binom_sf_mixture(n - (self.allele_tr - 1), n, p, w, smooth)) / norm_coef
        else:
            return self.binom_sf_mixture(x - 1, n, p, w, smooth)

    @staticmethod
    def odds_es(x, n, BAD, w):
        # w * np.log2(OR / BAD) + (1 - w) * np.log2(OR * BAD)
        return np.log2(x) - np.log2(n - x) + (1 - 2 * w) * np.log2(BAD) 
    
    @staticmethod
    def modify_w_binom(counts, n, p, ws):
        idx = (ws != 1) & (ws != 0)
        ws[idx] = expit(
            -(np.log(ws[idx]) - np.log(1 - ws[idx])
            + (n - 2 * counts[idx]) * (np.log(1 - p[idx]) - np.log(p[idx]))
            )
        )

        return ws


def calc_fdr(group_df):
    ind = group_df['min_pval'].notna()
    if group_df[ind].shape[0] > 0:  # check if any non-NA p-values exist
        corrected_pvalues = multipletests(
            group_df[ind]['min_pval'], alpha=0.05, method='fdr_bh')[1]
        group_df.loc[ind, 'FDR_sample'] = corrected_pvalues
    return group_df

def main(df, coverage_tr='auto', allele_tr=5, modify_w=False):
    df = df[df.eval(f'alt_counts >= {allele_tr} & ref_counts >= {allele_tr}')]
    # Remove already present columns
    df = df[[x for x in df.columns if x not in updated_columns]]
    # Check if empty
    result_columns = [*df.columns, 'coverage', *updated_columns]
    if df.empty:
        return pd.DataFrame([], columns=result_columns)

    imbalance_est = CalcImbalance(allele_tr=allele_tr, modify_w=modify_w)
    df['coverage'] = df.eval('ref_counts + alt_counts')
    if coverage_tr == 'auto':
        by_BAD_coverage_tr = {x: imbalance_est.calc_min_cover_by_BAD(x, pvalue_tr=0.05) for x in df['BAD'].unique()}
        df['is_tested'] = df['coverage'] >= df['BAD'].map(by_BAD_coverage_tr)
    else:    
        df['is_tested'] = df.eval(f'coverage >= {coverage_tr}')

    result = imbalance_est.calc_pval(
        df['coverage'].to_numpy(), 
        df['ref_counts'].to_numpy(),
        BADs=df['BAD'].to_numpy(),
        smooth=False
    )
    result = df.assign(**dict(zip(updated_columns, result)))[result_columns]
    result['min_pval'] = result[['pval_ref', 'pval_alt']].min(axis=1) * 2
    ind = (result['min_pval'] > 1) | (result['is_tested'])
    result['min_pval'] = np.where(ind, pd.NA, result['min_pval'])
    result['FDR_sample'] = pd.NA
    return result.groupby('sample_id').progress_apply(calc_fdr).reset_index()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('-a', type=int, help='Allelic reads threshold', default=5)
    parser.add_argument('--recalc-w', help='Specify to recalculate w',
        default=False, action="store_true")
    parser.add_argument('--ct', type=str, help="""Coverage threshold for individual variants to be considered tested.
                                        Expected to be "auto" or a positive integer""", default='auto')
    args = parser.parse_args()
    try:
        coverage_tr = int(args.ct) if args.ct != 'auto' else 'auto'
    except ValueError:
        print(f'Incorrect coverage threshold provided. {args.ct} not a positive integer or "auto"')
        raise
    input_df = pd.read_table(args.I, low_memory=False)
    modified_df = main(
        input_df,
        allele_tr=args.a,
        modify_w=args.recalc_w,
        coverage_tr=coverage_tr
    )
    modified_df.to_csv(args.O, sep='\t', index=None)