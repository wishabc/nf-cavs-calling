import pandas as pd
from scipy.stats import binom
from scipy.special import betainc
import argparse
import numpy as np


updated_columns = ['w', 'es', 'pval_ref', 'pval_alt', 'is_tested']

def calc_min_cover_by_BAD(BAD, es=1, pvalue_tr=0.05, allele_tr=5, cmax=1000):
    covs = np.arange(allele_tr * 2, cmax + 1)
    x = covs / (1 + np.power(2.0, -es) / BAD)
    mask = (allele_tr <= x) & (x <= covs - allele_tr)
    x = x[mask]
    covs = covs[mask]
    BADs = np.full(x.shape, BAD)
    _, _, sigs, _ = calc_pval(covs, x, BADs, allele_tr=allele_tr, modify_w=True, smooth=True)
    return np.min(covs[sigs < pvalue_tr])


def cdf(x, n, p, smooth):
    if smooth:
        return betainc(n - x, x + 1, 1 - p)
    else:
        return binom.cdf(x, n, p)


def sf(x, n, p, smooth):
    if smooth:
        return betainc(x + 1, n - x, p)        
    else:
        return binom.sf(x, n, p)


def binom_sf(x, n, p, w, smooth):
    return (1 - w) * sf(x, n, p, smooth) + w * sf(x, n, 1 - p, smooth)

def binom_cdf(x, n, p, w, smooth):
    return (1 - w) * cdf(x, n, p, smooth) + w * cdf(x, n, 1 - p, smooth)


def censored_binom_pvalue(x, n, p, w, allele_tr, smooth):
    if allele_tr > 0:
        norm_coef = binom_cdf(n - (allele_tr - 1), n, p, w, smooth) - binom_cdf(allele_tr - 1, n, p, w, smooth)
        return (binom_sf(x - 1, n, p, w, smooth) - binom_sf(n - (allele_tr - 1), n, p, w, smooth)) / norm_coef
    else:
        return binom_sf(x - 1, n, p, w, smooth)

def binom_es(x, n, p, w):
    return np.log2(
        x / (w * (n * (1 - p)) + (1 - w) * (n * p))
    )

def modify_w_binom(n, p, ws, counts):
    p1 = ws * binom.pmf(counts, n, 1 - p)
    p2 = (1 - ws) * binom.pmf(counts, n, p)
    return recalc_ws(ws, p1, p2)


# Both negbin and binom
def odds_es(x, n, BAD, w):
    return np.log2(x) - np.log2(n - x) + (1 - 2 * w) * np.log2(BAD) 
    
    # w * np.log2(OR / BAD) + (1 - w) * np.log2(OR * BAD)


def recalc_ws(ws, p1, p2):
    idx = (ws != 1) & (ws != 0)
    ws[idx] = p1[idx] / (p1[idx] + p2[idx])
    return ws

def calc_pval(n, k, BADs, allele_tr, modify_w, smooth=False):
    p = 1 / (BADs + 1)
    ws = np.full(k.shape, 0.5, dtype=np.float128)
    if modify_w:
        ws = modify_w_binom(n, p, ws, k)
    
    es = odds_es(k, n, BADs, ws)
    pval_ref = censored_binom_pvalue(k, n, p, ws, allele_tr, smooth=smooth)
    pval_alt = censored_binom_pvalue(n - k, n, p, 1 - ws, allele_tr, smooth=smooth)
    return [ws, es, pval_ref, pval_alt]
    
def calc_pval_for_df(df, allele_tr, modify_w):
    n = df.eval('alt_counts + ref_counts').to_numpy()
    ref_counts = df['ref_counts'].to_numpy()
    BADs = df['BAD'].to_numpy()

    return df.assign(
        **dict(zip(
            updated_columns,
            calc_pval(
                n, ref_counts,
                BADs=BADs,
                allele_tr=allele_tr,
                modify_w=modify_w,
                smooth=False
            )
        ))
    )


def main(df, coverage_tr='auto', allele_tr=5, modify_w=False):
    df = df[df.eval(f'alt_counts >= {allele_tr} & ref_counts >= {allele_tr}')]
    # Remove already present columns
    df = df[[x for x in df.columns if x not in updated_columns]]
    # Check if empty
    result_columns = [*df.columns, 'coverage', *updated_columns]
    if df.empty:
        return pd.DataFrame([], columns=result_columns)

    df['coverage'] = df.eval('ref_counts + alt_counts')
    if coverage_tr == 'auto':
        by_BAD_coverage_tr = {x: calc_min_cover_by_BAD(x) for x in df['BAD'].unique()}
        df['is_tested'] = df['coverage'] >= df['BAD'].apply(lambda x: by_BAD_coverage_tr[x])
    else:    
        df['is_tested'] = df.eval(f'coverage >= {coverage_tr}')

    return calc_pval_for_df(df, allele_tr=allele_tr, modify_w=modify_w)[result_columns]


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