import pandas as pd
from scipy.stats import binom
from scipy.special import betainc
import argparse
import numpy as np

from scipy.special import expit, logit
from aggregation import calc_fdr_pd, get_min_pval, logit_es, filter_pval_df

from tqdm import tqdm

tqdm.pandas()


updated_columns = ['coverage', 'w', 'es', 'logit_es', 'pval_ref', 'pval_alt', 'min_pval', 'FDR_sample']

class CalcImbalance:
    def __init__(self, allele_tr, modify_w):
        self.allele_tr = allele_tr
        self.modify_w = modify_w
    
    def calc_pval(self, n, k, BADs, smooth=False):
        p = 1 / (BADs + 1)
        ws = np.full(k.shape, 0.5, dtype=np.float64)
        if self.modify_w:
            ws = self.modify_w_binom(k, n, p, ws)
        
        fraction_es = self.odds_es(k, n, BADs, ws)

        pval_ref = self.censored_binom_pvalue(k, n, p, ws, smooth=smooth)
        pval_alt = self.censored_binom_pvalue(n - k, n, p, 1 - ws, smooth=smooth)
        return [ws, fraction_es, pval_ref, pval_alt]

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
        x, n, BAD = map(np.asarray, [x, n, BAD])

        b = np.log(BAD)
        delta = b * (n - 2 * x)
        p = x / n
        return expit(expit(-delta) * (logit(p) - b) + expit(delta) * (logit(p) + b))
    
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
    corrected_pvalues = calc_fdr_pd(group_df['min_pval'])
    group_df.loc[:, 'FDR_sample'] = corrected_pvalues
    return group_df


def main(df, coverage_tr=15, modify_w=False):
    df = df.drop(columns=updated_columns, errors='ignore')
    # Check if empty
    result_columns = [*df.columns, *updated_columns]
    if df.empty:
        return pd.DataFrame([], columns=result_columns)

    imbalance_est = CalcImbalance(allele_tr=0, modify_w=modify_w)
    df['coverage'] = df.eval('ref_counts + alt_counts')

    result = imbalance_est.calc_pval(
        df['coverage'].to_numpy(), 
        df['ref_counts'].to_numpy(),
        BADs=df['BAD'].to_numpy(),
        smooth=False
    )

    result = df.assign(
        **dict(zip(['w', 'es', 'pval_ref', 'pval_alt'], result))
    )
    result['logit_es'] = logit_es(result['es'])
    result = result.groupby('sample_id').progress_apply(filter_pval_df, max_cover_tr=coverage_tr)
    result['min_pval'] = get_min_pval(
        result,
        columns=['pval_ref', 'pval_alt']
    )
    result['FDR_sample'] = result.groupby('sample_id')['min_pval'].transform(calc_fdr_pd)
    return result.reset_index(drop=True)[result_columns]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('--recalc-w', help='Specify to recalculate w',
        default=False, action="store_true")
    parser.add_argument('--coverage_threhold', type=int, help="""Coverage threshold for variants to calculate per-sample q-values.
                            Expected to a positive integer""", default=15)
    args = parser.parse_args()
    input_df = pd.read_table(args.I, low_memory=False)
    modified_df = main(
        input_df,
        modify_w=args.recalc_w,
        coverage_tr=args.coverage_threhold
    )
    modified_df.to_csv(args.O, sep='\t', index=None)