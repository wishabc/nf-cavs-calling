import pandas as pd
from scipy.stats import binom
from scipy.special import betainc
import argparse
import numpy as np


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
    return ws, es, pval_ref, pval_alt

    
def calc_pval_for_df(df, allele_tr, modify_w):
    n = df.eval('alt_counts + ref_counts').to_numpy()
    ref_counts = df['ref_counts'].to_numpy()
    BADs = df['BAD'].to_numpy()

    ws, es, pval_ref, pval_alt = calc_pval(
        n, ref_counts, BADs,
        allele_tr=allele_tr,
        modify_w=modify_w,
        smooth=False
    )
    df['w'] = ws
    df['es'] = es
    df['pval_ref'] = pval_ref
    df['pval_alt'] = pval_alt
    return df


def main(df, allele_tr=5, modify_w=False):
    df = df[df.eval(f'alt_counts >= {allele_tr} & ref_counts >= {allele_tr}')]
    # Remove already present columns
    df = df[[x for x in df.columns if x not in ['w', 'es', 'pval_ref', 'pval_alt']]]
    # Check if empty
    if df.empty:
        cols = ['w', 'es', 'pval_ref', 'pval_alt']
        for col in cols:  
            df[col] = np.nan
    else:
        df = calc_pval_for_df(df, allele_tr=allele_tr, modify_w=modify_w)
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('-a', type=int, help='Allelic reads threshold', default=5)
    parser.add_argument('--recalc-w', help='Specify to recalculate w',
        default=False, action="store_true")
    args = parser.parse_args()

    input_df = pd.read_table(args.I)
    modified_df = main(
        input_df,
        allele_tr=args.a,
        modify_w=args.recalc_w
    )
    modified_df.to_csv(args.O, sep='\t', index=None)