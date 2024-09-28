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

def calc_pval(self, n, k, BADs, smooth=False):
    p = 1 / (BADs + 1)
    ws = np.full(k.shape, 0.5, dtype=np.float64)
    if self.modify_w:
        ws = self.modify_w_binom(k, n, p, ws)
    
    fraction_es = self.odds_es(k, n, BADs, ws)

    pval_ref = self.censored_binom_pvalue(k, n, p, ws, smooth=smooth)
    pval_alt = self.censored_binom_pvalue(n - k, n, p, 1 - ws, smooth=smooth)
    return [ws, fraction_es, pval_ref, pval_alt]

def odds_es(x, n, BAD, w):
    # w * np.log2(OR / BAD) + (1 - w) * np.log2(OR * BAD)
    x, n, BAD = map(np.asarray, [x, n, BAD])

    log_bad = np.log(BAD)
    delta = log_bad * (n - 2 * x)
    p = x / n
    return expit(expit(-delta) * (logit(p) - log_bad) + expit(delta) * (logit(p) + log_bad))
    


def calc_fdr(group_df):
    corrected_pvalues = calc_fdr_pd(group_df['min_pval'])
    group_df.loc[:, 'FDR_sample'] = corrected_pvalues
    return group_df


def main(df, coverage_tr=15):
    df = df.drop(columns=updated_columns, errors='ignore')
    # Check if empty
    result_columns = [*df.columns, *updated_columns]
    if df.empty:
        return pd.DataFrame([], columns=result_columns)

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
    result = result.groupby('sample_id').progress_apply(filter_pval_df, max_cover_tr=coverage_tr).reset_index(drop=True)
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
    parser.add_argument('--coverage_threhold', type=int, help="""Coverage threshold for variants to calculate per-sample q-values.
                        Expected to be a positive integer""", default=15)
    args = parser.parse_args()
    input_df = pd.read_table(args.I, low_memory=False)
    modified_df = main(
        input_df,
        coverage_tr=args.coverage_threhold
    )
    modified_df.to_csv(args.O, sep='\t', index=None)