from re import I
import sys
import os
import pandas as pd
from scipy.stats import binom, nbinom
from .helpers import alleles
import argparse


result_columns = ['#chr', 'start', 'end', 'ID', 'ref', 'alt', 'ref_counts', 'alt_counts', 'BAD', 'pval_ref', 'pval_alt']


def nbinom_sf(x, r, p, w):
    return w * nbinom.sf(x, r, p) + (1 - w) * nbinom.sf(x, r, 1 - p)


def censored_nbinom_pvalue(x, r, p, w, allele_tr):
    norm_coef = 0
    if allele_tr > 0:
        norm_coef = nbinom_sf(allele_tr - 1, r, p, w)
    else:
        norm_coef = 1
    return nbinom_sf(x - 1, r, p, w) / norm_coef

def binom_sf(x, n, p):
    return 0.5 * (binom.sf(x, n, p) + binom.sf(x, n, 1 - p))

def binom_cdf(x, n, p):
    return 0.5 * (binom.cdf(x, n, p) + binom.cdf(x, n, 1 - p))

def censored_binom_pmf(n, p, allele_tr):
    dist1 = binom(n, p)
    dist2 = binom(n, 1 - p)
    norm = (dist1.cdf(n - allele_tr) + dist2.cdf(n - allele_tr)
    - dist1.cdf(allele_tr - 1) - dist2.cdf(allele_tr - 1))
    return [
        (dist1.pmf(i) + dist2.pmf(i)) / norm if (i >= allele_tr and i <= n - allele_tr) else 0
        for i in range(allele_tr, n - allele_tr + 1)
    ]


def censored_binom_pvalue(x, n, p, allele_tr):
    norm_coef = 0
    if allele_tr > 0:
        norm_coef = binom_cdf(n - (allele_tr - 1), n, p) - binom_cdf(allele_tr - 1, n, p)
    else:
        norm_coef = 1
    return (binom_sf(x - 1, n, p) - binom_sf(n - (allele_tr - 1), n, p)) / norm_coef


def calc_pval_for_indiv(input_filename, output_filename, stats_file=None, mode='binom', allele_tr=5):
    df = pd.read_table(input_filename)
    if stats_file is not None:
        stats_df = pd.read_table(stats_file)
    else:
        stats_df = None
    df = calc_pval_for_df(df, stats_df, mode, allele_tr)
    df = df[result_columns]
    df.to_csv(output_filename, sep='\t', index=None)
        
def calc_pval_for_df(df, nb_params, mode='binom', allele_tr=5):
    p = df.eval('1 / (BAD + 1)').to_numpy()
    n = df.eval('alt_counts + ref_counts').to_numpy()
    for allele in alleles:
        pval_field = f'pval_{allele}'
        if mode == 'binom':
            df[pval_field] = censored_binom_pvalue(df[f'{allele}_counts'].to_numpy(), n, p, allele_tr)
        elif mode == 'negbin':
            if nb_params is None:
                raise ValueError('NB params are required for p-value calculations')
            merged = df.reset_index(drop=True).reset_index()\
            .merge(nb_params[nb_params.eval(f'allele == "{allele}"')],
             left_on=[f'{alleles[allele]}_counts', 'BAD'],
             right_on=['fix_c', 'BAD'], sort=False, how='left')\
            .sort_values('index')
            # if conservative:
            #     merged['w'] = 1
            #     merged['r'] = merged['fix_c']
            # else:
            merged.loc[merged.eval('r == 0'), 'w'] = 1
            merged.loc[merged.eval('r == 0'), 'r'] = merged.loc[merged.eval('r == 0'), 'fix_c']
            rs = merged['r'].to_numpy()
            ws = merged['w'].to_numpy()
            df[pval_field] = censored_nbinom_pvalue(df[f'{allele}_counts'].to_numpy(), rs, p, ws, allele_tr)
        else:
            raise AssertionError
    return df


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('-s', help='Strategy for p-value calculation. One of "negbin", "binom"')
    parser.add_argument('--stats-file',
        help='Path to file with fitted negative binomial distribution. Not required for "binom" strategy')
    parser.add_argument('-a', type=int, help='Allelic reads threshold', default=5)
    args = parser.parse_args()
    calc_pval_for_indiv(
        input_filename=args['-I'],
        output_filename=args['-O'],
        mode=args['-s'],
        stats_file=args['--stats_file'],
        allele_tr=args['-a']
    )