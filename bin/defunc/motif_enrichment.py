#!/bin/which python3

import sys
import pandas as pd
import numpy as np
import scipy.stats as stats

# Helper functions

_complement = {"A": "T", "C": "G", "G": "C", "T": "A"}

def complement(x):
    return _complement[x]

def reverse_complement(x):    
    bases = list(x) 
    bases = reversed([_complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases

def flip_by_strand(x):
    if x['strand'] == '-':
        x['ref'] = complement(x['ref'])
        x['alt'] = complement(x['alt'])
    return x

def prefered_allele(x):
    x["prefered_allele"] = x["ref"] if x["ard"] >= 0.5 else x["alt"]
    return x


# Load variant imbalance file
variants_df = pd.read_table(sys.argv[1], header=0)
variants_df.set_index('variant_id', inplace=True)
variants_df = variants_df[~np.isnan(variants_df['adj_p'])]

# Load motifs dataframe
motifs_df = pd.read_table(sys.argv[2], header=None)
motifs_df.columns = ['chrom', 'start', 'end', 'variant_id', 'dbsnp', 'ref', 'alt', 'aa', 'motif', 'offset', 'within', 'strand', 'ref_score', 'alt_score', 'seq']
motifs_df.set_index('variant_id', inplace=True)

# Flip data with motif is on '-' strand
motifs_df = motifs_df.apply(flip_by_strand, axis = 1)

# Add imbalance data
motifs_df = motifs_df.merge(variants_df[['ard','adj_p']].reset_index(), on='variant_id')

# Compute preferred allele
motifs_df = motifs_df.apply(prefered_allele, axis = 1)

# Params
flank_width = 20
n_shuffles = 1000

#Enrichment code

imbalanced = (motifs_df['adj_p'] <= 0.01) & ((motifs_df['ard'] > 0.65))# | (df['ard'] < 0.3))

w = motifs_df['offset'].max()+1-flank_width

bins = np.arange(motifs_df['offset'].min(), motifs_df['offset'].max()+2) # add 2 to length: one for '0' and one for last element

n_all = np.histogram(motifs_df['offset'], bins=bins)[0]
n_imbalanced = np.histogram(motifs_df['offset'][imbalanced], bins=bins)[0]
n_not_imbalanced = n_all - n_imbalanced

log_odds = np.log2( (n_imbalanced[flank_width:-flank_width].sum() / n_imbalanced.sum()) / (n_not_imbalanced[flank_width:-flank_width].sum() / n_not_imbalanced.sum()) )
log_odds_per_nt = np.log2( (n_imbalanced / n_imbalanced.sum()) / (n_not_imbalanced / n_not_imbalanced.sum()) )

perm = np.zeros(n_shuffles)
perm_per_nt = np.zeros((n_shuffles, len(bins)-1))

for i in range(n_shuffles):
    n_exp_imbalanced = np.histogram(motifs_df['offset'][np.random.permutation(imbalanced)], bins=bins)[0] + 1
    n_exp_not_imbalanced = n_all - n_exp_imbalanced +1 

    perm[i] = np.log2( (n_exp_imbalanced[flank_width:-flank_width].sum() / n_exp_imbalanced.sum()) / (n_exp_not_imbalanced[flank_width:-flank_width].sum() / n_exp_not_imbalanced.sum()) )
    perm_per_nt[i,:] = np.log2( (n_exp_imbalanced / np.sum(n_exp_imbalanced)) / (n_exp_not_imbalanced / np.sum(n_exp_not_imbalanced)))

pval = -stats.norm.logsf(log_odds, loc=np.nanmean(perm, axis=0), scale=np.nanstd(perm, axis=0))/np.log(10)
pvals_per_nt = -stats.norm.logsf(log_odds_per_nt, loc=np.nanmean(perm_per_nt, axis=0), scale=np.nanstd(perm_per_nt, axis=0))/np.log(10)

motif, cluster = sys.argv[3].split(':')

fields = [
    motif,
    cluster,
    log_odds,
    pval,
    np.nansum(n_all[flank_width:-flank_width]),
    np.nansum(n_imbalanced[flank_width:-flank_width]),
    np.nanmedian(n_all[flank_width:-flank_width]),
    np.nansum(n_imbalanced[flank_width:-flank_width]>=7)
]

print('\t'.join(map(str, fields)))