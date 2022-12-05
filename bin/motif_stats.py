import pandas as pd
import scipy as sp

import numpy as np
import statsmodels.api as sm
from tqdm import tqdm
import sys
import os


tqdm.pandas()


_complement = {"A": "T", "C": "G", "G": "C", "T": "A"}

def complement(x):
    return _complement[x]


def set_index(df):
    df.index = df['#chr'] + ':' + df['start'].astype(str) + ':' + df['alt']
    return df


def flip_by_strand(x):
    if x['strand'] == '-':
        x['ref'] = complement(x['ref'])
        x['alt'] = complement(x['alt'])
        #x['aa'] =  complement(x['aa'])
    return x

def prefered_allele(x, es_fld='aggregated_es_weighted_ref'):
    x["prefered_allele"] = x["ref"] if x[es_fld] >= 0 else x["alt"]
    return x

def logtr(arr):
    return 1 - 1 / (1 + np.power(2, np.array(arr)))

def get_concordant(x, y, x0=0, x_mar=0):
    return ((y * (x - x0) > 0) & (np.abs(x - x0) >= x_mar)).sum() / ((np.abs(x - x0) >= x_mar) & (y != 0)).sum()

def variant_enrichment(variant_by_motif, n_shuffles=1000):
    
    df = variant_by_motif.copy()
    imbalanced = (df['min_fdr'] <= 0.01) #& ((df['ard'] > 0.6))# | (df['ard'] < 0.3))

    flank_width = 20
    w = df['offset'].max()+1-flank_width
    motif_positions = np.arange(w)
    
    bins = np.arange(df['offset'].min(), df['offset'].max()+2) # add 2 to length: one for '0' and one for last element

    n_all = np.histogram(df['offset'], bins=bins)[0]
    n_imbalanced = np.histogram(df['offset'][imbalanced], bins=bins)[0]
    n_not_imbalanced = n_all - n_imbalanced
    
    log_odds = np.log2( (n_imbalanced[flank_width:-flank_width].sum() / n_imbalanced.sum()) / (n_not_imbalanced[flank_width:-flank_width].sum() / n_not_imbalanced.sum()) )
    log_odds_per_nt = np.log2( (n_imbalanced / n_imbalanced.sum()) / (n_not_imbalanced / n_not_imbalanced.sum()) )

    perm = np.zeros(n_shuffles)
    perm_per_nt = np.zeros((n_shuffles, len(bins)-1))

    for i in range(n_shuffles):
        n_exp_imbalanced = np.histogram(df['offset'][np.random.permutation(imbalanced)], bins=bins)[0] + 1
        n_exp_not_imbalanced = n_all - n_exp_imbalanced +1 

        perm[i] = np.log2( (n_exp_imbalanced[flank_width:-flank_width].sum() / n_exp_imbalanced.sum()) / (n_exp_not_imbalanced[flank_width:-flank_width].sum() / n_exp_not_imbalanced.sum()) )
        if np.sum(n_exp_imbalanced) == 0 or np.sum(n_exp_not_imbalanced) == 0:
            perm_per_nt[i,:] = np.nan
        else:
            perm_per_nt[i,:] = np.log2( (n_exp_imbalanced / np.sum(n_exp_imbalanced)) / (n_exp_not_imbalanced / np.sum(n_exp_not_imbalanced)))
  
    pval = -sp.stats.norm.logsf(log_odds, loc=np.nanmean(perm, axis=0), scale=np.nanstd(perm, axis=0))/np.log(10)
    pvals_per_nt = -sp.stats.norm.logsf(log_odds_per_nt, loc=np.nanmean(perm_per_nt, axis=0), scale=np.nanstd(perm_per_nt, axis=0))/np.log(10)
    
    return df, imbalanced, log_odds, pval, log_odds_per_nt, pvals_per_nt, n_imbalanced[flank_width:-flank_width].sum(), n_all[flank_width:-flank_width].sum()


def get_annotations(motif_id, variants, motif_counts):
    df = pd.read_table(motif_counts, header=None,
                       names=['#chr', 'start', 'end', 'rsID',
                        'ref', 'alt', 'motif', 'offset',
                        'within', 'strand', 'ref_score', 'alt_score', 'seq'])
    v = set_index(df)
    v = v.apply(flip_by_strand, axis = 1)
    es_flds = ['aggregated_es_weighted_ref', 'aggregated_es_ref']
    snvs = v.join(variants[es_flds + ['min_fdr']])
    snvs = snvs.apply(lambda x: prefered_allele(x, 'aggregated_es_weighted_ref'), axis=1)
    var, imbl, log_odds, pval, _, pvals_per_nt, n_imb, n_all = variant_enrichment(snvs)

    df_ct_all = pd.crosstab(var['offset'], var['ref'])
    df_ct_imbalanced = pd.crosstab(var.loc[imbl]['offset'], var.loc[imbl]['prefered_allele'])
   
    inside = snvs['within'] == 1
    imb = snvs['min_fdr']<=0.01
    
    snvs_in = snvs[inside & imb]

    pred = 'log_es'
    es_fld = 'aggregated_es_weighted_ref'

    if pred == 'log_es':
        lim = (-3, 3)
        nbins = lim[1] * 8 + 1
        fit_lim = (-2, 2)
        predictor_array = snvs_in[es_fld]
        x_0 = 0
    elif pred == 'raw_es':
        lim = (0, 1)
        nbins = 25
        predictor_array = logtr(snvs_in[es_fld])
        fit_lim = (0.2, 0.8)
        x_0 = 0.5
    elif pred == 'signed_fdr':
        lim = (-25, 25)
        nbins = 25
        predictor_array = -np.log10(snvs_in['min_fdr']) * np.sign(snvs_in['aggregated_es_weighted_ref'])
        fit_lim = (-15, 15)
        x_0 = 0

    snvs_in['ddg'] = snvs_in['ref_score'] - snvs_in['alt_score']

    bins=np.linspace(*lim, nbins)

    xdiff = (bins[1]-bins[0])/2.0
    x = bins[:-1]+xdiff

    snvs_in["bin"] = pd.cut(predictor_array, bins)
    
    
    X = predictor_array
    if n_imb > 0:
            
        Y = snvs_in['ddg']

        X = sm.add_constant(X)


        model = sm.OLS(Y, X).fit()

        outliers = model.get_influence().cooks_distance[0] > (4/X.shape[0])

        model = sm.OLS(Y[~outliers], X[~outliers]).fit()

        xx = np.array(fit_lim)
        yy = xx*model.params.iloc[1] + model.params.iloc[0]
        
        r2 = model.rsquared
        conc = get_concordant(predictor_array, snvs_in["ddg"], x0=x_0)
        log_conc = np.log2(conc) - np.log2(1 - conc)
        
        log_ref_bias = np.log2((predictor_array > 0).sum()) - np.log2((predictor_array <0).sum())
    else:
         r2, log_conc, log_ref_bias = np.nan, np.nan, np.nan

    return pd.DataFrame.from_records(
        [[motif_id, log_odds, pval, n_imb, n_all, r2, log_conc, log_ref_bias]],
        columns=['motif_id', 'log_odds', 'pval','n_imb', 'n_all', 
        'r2', 'log_conc', 'log_ref_bias'])
        
    
    
def main(motif_id, all_variants, motif_counts, out_path):
    variants = pd.read_table(all_variants, header=0)
    variants['min_fdr'] = variants[['fdrp_bh_ref', 'fdrp_bh_alt']].min(axis=1)
    variants = set_index(variants)

    motif_counts = get_annotations(motif_id=motif_id, motif_counts=motif_counts, variants=variants)

    motif_counts.to_csv(out_path, sep='\t', index=False)

if __name__ == '__main__':
    main(*sys.argv[1:])