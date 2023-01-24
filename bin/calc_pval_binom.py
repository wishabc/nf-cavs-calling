import pandas as pd
from scipy.stats import binom, nbinom
from helpers import alleles
import argparse
import numpy as np


# nbinom_dist
def nbinom_sf(x, r, p, w):
    return w * nbinom.sf(x, r, p) + (1 - w) * nbinom.sf(x, r, 1 - p)


def nbinom_pmf(x, r, p, w):
    return w * nbinom.pmf(x, r, p) + (1 - w) * nbinom.pmf(x, r, 1 - p)

def map_func_to_array(arr, func):
    return arr * func(arr)

def censored_nbinom_expectation(x, r, p, w, allele_tr):
    #vector_map = np.vectorize(map_func_to_array)
    return (r * ((1 - p) / p * w + (1 - w) * p / (1 - p)) - map_func_to_array(np.tile(np.arange(5), (x.shape[0], 1)).transpose(), lambda y: nbinom_pmf(y, r, p, w)).sum(axis=0)) / nbinom_sf(allele_tr - 1, r, p, w)
    
def censored_nbinom_es(x, r, p, w, allele_tr):
    return np.log2(x / censored_nbinom_expectation(x, r, p, w, allele_tr))

def censored_nbinom_pvalue(x, r, p, w, allele_tr):
    norm_coef = 0
    if allele_tr > 0:
        norm_coef = nbinom_sf(allele_tr - 1, r, p, w)
    else:
        norm_coef = 1

    return nbinom_sf(x - 1, r, p, w) / norm_coef

def modify_w_nbinom(rs, p, ws, counts):
    p1 = ws * nbinom.pmf(counts, rs, p)
    p2 = (1 - ws) * nbinom.pmf(counts, rs, 1 - p)
    return recalc_ws(ws, p1, p2)

# binom dist
def binom_sf(x, n, p, w):
    return (1 - w) * binom.sf(x, n, p) + w * binom.sf(x, n, 1 - p)

def binom_cdf(x, n, p, w):
    return (1 - w) * binom.cdf(x, n, p) + w * binom.cdf(x, n, 1 - p)

def binom_pmf(x, n, p, w):
    return (1 - w) * binom.pmf(x, n, p) + w * binom.pmf(x, n, 1 - p)

def censored_binom_pmf(n, p, allele_tr):
    dist1 = binom(n, p)
    dist2 = binom(n, 1 - p)
    norm = (dist1.cdf(n - allele_tr) + dist2.cdf(n - allele_tr)
    - dist1.cdf(allele_tr - 1) - dist2.cdf(allele_tr - 1))
    return [
        (dist1.pmf(i) + dist2.pmf(i)) / norm if (i >= allele_tr and i <= n - allele_tr) else 0
        for i in range(allele_tr, n - allele_tr + 1)
    ]

def censored_binom_pvalue(x, n, p, w, allele_tr):
    norm_coef = 0
    if allele_tr > 0:
        norm_coef = binom_cdf(n - (allele_tr - 1), n, p, w) - binom_cdf(allele_tr - 1, n, p, w)
    else:
        norm_coef = 1
    return (binom_sf(x - 1, n, p, w) - binom_sf(n - (allele_tr - 1), n, p, w)) / norm_coef


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



def calc_pval_for_indiv(input_filename, output_filename, allele_tr=5, modify_w=False, es_method='exp'):
    df = pd.read_table(input_filename)
    df = calc_pval_for_df(df, allele_tr,modify_w=modify_w, es_method=es_method)
    df.to_csv(output_filename, sep='\t', index=None)


def calc_pval_for_df(df, allele_tr, modify_w, es_method):
    if df.empty:
        cols = ['w', 'es', 'pval_ref', 'pval_alt']
        for col in cols:  
            df[col] = np.nan
        return df
    p = df.eval('1 / (BAD + 1)').to_numpy()
    n = df.eval('alt_counts + ref_counts').to_numpy()
    ref_counts = df['ref_counts'].to_numpy()
    BADS = df['BAD'].to_numpy()

    ws = np.full(ref_counts.shape, 0.5, dtype=np.float128)
    if modify_w:
        ws = modify_w_binom(n, p, ws, ref_counts)
    df['w'] = ws
    if es_method == 'exp':
        es_list = binom_es(ref_counts, n, p, ws)
    elif es_method == 'odds':
        es_list = odds_es(ref_counts, n, BADS, ws)
    elif es_method == 'cons':
        es_list = odds_es(ref_counts, n, BADS, 1)
    df['es'] = es_list

    df['pval_ref'] = censored_binom_pvalue(ref_counts, n, p, ws, allele_tr)
    df['pval_alt'] = censored_binom_pvalue(n - ref_counts, n, p, 1 - ws, allele_tr)
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('-a', type=int, help='Allelic reads threshold', default=5)
    parser.add_argument('--es-method',
        help='Method to calculate effect size. One of "exp", "odds", "cons"',
        default='odds')
    parser.add_argument('--recalc-w', help='Specify to recalculate w',
        default=False, action="store_true")
    args = parser.parse_args()
    calc_pval_for_indiv(
        input_filename=args.I,
        output_filename=args.O,
        allele_tr=args.a,
        es_method=args.es_method,
        modify_w=args.recalc_w
    )