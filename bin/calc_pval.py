import pandas as pd
from scipy.stats import binom, nbinom
from helpers import alleles, get_field_by_ftype
import argparse
import numpy as np

result_ftypes = ('es', 'pval')
result_columns = ['#chr', 'start', 'end', 'ID', 'ref', 'alt', 'ref_counts', 'alt_counts', 'sample_id', 'BAD'] + [get_field_by_ftype(allele, ftype) for allele in alleles for ftype in result_ftypes]


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
def odds_es(x, n, p, w):
    return np.log2(
        x / ((n - x) * ((1 - p) / p * w + (1 - w) * p / (1 - p)))
    )

def recalc_ws(ws, p1, p2):
    idx = (ws != 1) & (ws != 0)
    ws[idx] = p1[idx] / (p1[idx] + p2[idx])
    return ws



def calc_pval_for_indiv(input_filename, output_filename, stats_file=None, mode='binom', allele_tr=5, modify_w=False, es_method='exp'):
    df = pd.read_table(input_filename)
    if mode == 'negbin' and stats_file is not None:
        stats_df = pd.read_table(stats_file)
    else:
        stats_df = None
    df = calc_pval_for_df(df, stats_df, mode, allele_tr,modify_w=modify_w, es_method=es_method)
    df[result_columns].to_csv(output_filename, sep='\t', index=None)


def calc_pval_for_df(df, nb_params, mode, allele_tr, modify_w, es_method):
    if df.empty:
        for allele in alleles:
            for ftype in result_ftypes:
                df[get_field_by_ftype(allele, ftype=ftype)] = None
        return df
    p = df.eval('1 / (BAD + 1)').to_numpy()
    n = df.eval('alt_counts + ref_counts').to_numpy()
    for allele in alleles:
        counts = df[f'{allele}_counts'].to_numpy()
        if mode == 'binom':
            ws = np.full(counts.shape, 1, dtype=np.float128)
            if modify_w:
                ws = modify_w_binom(n, p, ws, counts)
            df[f'w_{allele}'] = ws
            df[get_field_by_ftype(allele)] = censored_binom_pvalue(counts, n, p, ws, allele_tr)

            if es_method == 'exp':
                es_list = binom_es(counts, n, p, ws)
            elif es_method == 'odds':
                es_list = odds_es(counts, n, p, ws)
            elif es_method == 'cons':
                es_list = odds_es(counts, n, p, 1)
            df[get_field_by_ftype(allele, 'es')] = es_list
        elif mode == 'negbin':
            if nb_params is None:
                raise ValueError('NB params are required for p-value calculations')
            merged = df.reset_index(drop=True).reset_index()\
            .merge(nb_params[nb_params.eval(f'allele == "{allele}"')],
             left_on=[f'{alleles[allele]}_counts', 'BAD'],
             right_on=['fix_c', 'BAD'], sort=False, how='left', validate='m:1')\
            .sort_values('index')
            merged.loc[merged.eval('gof > 0.05'), 'r'] = 0
            merged.loc[merged.eval('r == 0'), 'w'] = 1
            merged.loc[merged.eval('r == 0'), 'r'] = merged.loc[merged.eval('r == 0'), 'fix_c']
            rs = merged['r'].to_numpy()
            ws = merged['w'].to_numpy()
            if modify_w:
                ws = modify_w_nbinom(rs, p, ws, counts)
            df[get_field_by_ftype(allele)] = censored_nbinom_pvalue(counts, rs, p, ws, allele_tr)
            if es_method == 'exp':
                es_list = censored_nbinom_es(counts, rs, p, ws, allele_tr)
            elif es_method == 'odds':
                es_list = odds_es(counts, n, p, ws)
            elif es_method == 'cons':
                es_list = odds_es(counts, n, p, 1)
            df[get_field_by_ftype(allele, 'es')] = es_list
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('-s', help='Strategy for p-value calculation. One of "negbin", "binom"')
    parser.add_argument('--stats-file',
        help='Path to file with fitted negative binomial distribution. Not required for "binom" strategy',
        default='./')
    parser.add_argument('-a', type=int, help='Allelic reads threshold', default=5)
    parser.add_argument('--es-method',
        help='Method to calculate effect size. One of "exp", "odds", "cons"',
        default='exp')
    parser.add_argument('--recalc-w', help='Specify to recalculate w',
        default=False, action="store_true")
    args = parser.parse_args()
    calc_pval_for_indiv(
        input_filename=args.I,
        output_filename=args.O,
        mode=args.s,
        stats_file=args.stats_file,
        allele_tr=args.a,
        es_method=args.es_method,
        modify_w=args.recalc_w
    )