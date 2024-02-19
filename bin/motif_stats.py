#!/bin/which python3
import argparse
import pandas as pd
import numpy as np
import scipy.stats as stats
from tqdm import tqdm
tqdm.pandas()
from aggregation import starting_columns

_complement = {"A": "T", "C": "G", "G": "C", "T": "A"}


class NoDataException(Exception):
    pass


class MotifEnrichment:
    columns = [
        "log_odds",
        "log10_pval",
        "total",
        "imbalanced",
        "total_inside",
        "imbalanced_inside",
        "imbalanced_inside_median",
        "n_imbalanced_more_7",
    ]
    def __init__(self, flank_width=20, n_shuffles=1000):
        self.flank_width = flank_width
        self.n_shuffles = n_shuffles

    
    def calc_log_odds(self, imbalanced, not_imbalanced):
        return np.log2(imbalanced[self.flank_width:-self.flank_width].sum()) \
            - np.log2(imbalanced.sum()) \
            - np.log2(not_imbalanced[self.flank_width:-self.flank_width].sum()) \
            + np.log2(not_imbalanced.sum())


    def calc_enrichment(self, df, imbalanced):
        bins = np.arange(df['offset'].min(), df['offset'].max() + 2) # add 2 to length: one for '0' and one for last element

        n_all = np.histogram(df['offset'], bins=bins)[0]
        n_imbalanced = np.histogram(df['offset'][imbalanced], bins=bins)[0] + 1
        n_not_imbalanced = n_all - n_imbalanced + 1
        
        log_odds = self.calc_log_odds(n_imbalanced, n_not_imbalanced)
        log_odds_per_nt = np.log2( (n_imbalanced / n_imbalanced.sum()) / (n_not_imbalanced / n_not_imbalanced.sum()) )

        perm = np.zeros(self.n_shuffles)
        perm_per_nt = np.zeros((self.n_shuffles, len(bins)-1))
        
        for i in range(self.n_shuffles):
            rng = np.random.default_rng(seed=i * self.n_shuffles * 10)
            n_exp_imbalanced = np.histogram(df['offset'][rng.permutation(imbalanced)], bins=bins)[0] + 1
            n_exp_not_imbalanced = n_all - n_exp_imbalanced + 1

            perm[i] = self.calc_log_odds(n_exp_imbalanced, n_exp_not_imbalanced)
            perm_per_nt[i,:] = np.log2( (n_exp_imbalanced / np.sum(n_exp_imbalanced)) / (n_exp_not_imbalanced / np.sum(n_exp_not_imbalanced)))
    

        pval = -stats.norm.logsf(
            log_odds,
            loc=np.nanmean(perm, axis=0),
            scale=np.nanstd(perm, axis=0)
        ) / np.log(10)

        pvals_per_nt = -stats.norm.logsf(
            log_odds_per_nt, 
            loc=np.nanmean(perm_per_nt, axis=0),
            scale=np.nanstd(perm_per_nt, axis=0)
        ) / np.log(10)
        
        n_imbalanced_inside = n_imbalanced[self.flank_width:-self.flank_width]
        return pd.Series([
            log_odds,
            pval,
            np.nansum(n_all),
            np.nansum(n_imbalanced),
            np.nansum(n_all[self.flank_width:-self.flank_width]),
            np.nansum(n_imbalanced_inside),
            np.nanmedian(n_imbalanced_inside),
            np.nansum(n_imbalanced_inside >= 7)
        ], index=self.columns), (pvals_per_nt, log_odds_per_nt), (perm, perm_per_nt)


    # wrapper of calc_enrichment to handle no data
    def get_group_stats(self, group_df):
        imbalanced_index = group_df['imbalanced']
        try:
            if imbalanced_index.sum() == 0:
                raise NoDataException()

            data, _, _ = self.calc_enrichment(group_df, imbalanced_index)
        except NoDataException:
            data = pd.Series(np.full(len(self.columns), pd.NA), index=self.columns)
        return data


    def get_motif_stats(self, data_df):
        return data_df.groupby(['motif_id', 'group_id']).progress_apply(self.get_group_stats).reset_index()


def preprocess_dfs(variants_df, motifs_df):
    # Compute preferred allele
    data_df = variants_df[[*starting_columns, 'logit_es_combined', 'group_id', 'min_fdr']].merge(
        motifs_df,
        on=starting_columns
    )
    for allele in ('ref', 'alt'):
        data_df[f'{allele}_by_motif'] = np.where(data_df['strand'] == '-', data_df[allele].map(_complement), data_df[allele])

    data_df["prefered_allele"] = np.where(
        data_df['logit_es_combined'] >= 0,
        data_df["ref_by_motif"],
        data_df["alt_by_motif"])
    data_df['ddg'] = data_df.eval('ref_score - alt_score')
    return data_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converts matrix in txt format to npy')
    parser.add_argument('variants', help='Path to file with CAVs significance')
    parser.add_argument('motifs', help='Path to file with motif hits at variants')
    parser.add_argument('outpath', help='Path to output file')
    parser.add_argument('--flank_width', help='Width of flanks', type=int, default=20)
    parser.add_argument('--fdr', help='FDR threshold for CAVs', type=float, default=0.1)
    args = parser.parse_args()

    print('Reading variants df')
    variants_df = pd.read_table(args.variants)
    variants_df = variants_df[variants_df['min_fdr'].notna()]
    
    print('Reading motifs df')
    motifs_df = pd.read_table(args.motifs).rename(columns={'motif': 'motif_id'})
    print('Adding fields')
    data_df = preprocess_dfs(variants_df, motifs_df)
    data_df['imbalanced'] = data_df['min_fdr'] <= args.fdr

    me = MotifEnrichment(flank_width=args.flank_width, n_shuffles=1000)
    me.get_motif_stats(data_df).to_csv(args.outpath, index=False, sep='\t')
