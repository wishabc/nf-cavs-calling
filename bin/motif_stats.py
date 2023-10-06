#!/bin/which python3

import argparse
import pandas as pd
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
from tqdm import tqdm
# import dask.dataframe as dd
tqdm.pandas()

# Params
_complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
def complement(x):
    return np.vectorize(_complement.get)(x)


class NoDataException(Exception):
    pass

class MotifEnrichment:
    result_columns =  [
        'motifs_name',
        'group_id',
        "log_odds",
        "pval",
        "total_inside",
        "imbalanced_inside",
        "imbalanced_inside_median",
        "n_imbalanced_more_7",
        'r2',
        'concordance',
        'ref_bias'
    ]

    def __init__(self, variants_df_path, counts_df_path, flank_width=20, fdr_tr=0.05):
        self.flank_width = flank_width
        self.fdr_tr = fdr_tr
        print('Reading variants df')
        variants_df = pd.read_table(variants_df_path).set_index(['#chr', 'start', 'end', 'ref', 'alt'])
        
        print('Reading motifs df')
        motifs_df = pd.read_table(counts_df_path).set_index(['#chr', 'start', 'end', 'ref', 'alt'])
        print('Adding fields')
        for key in ('ref', 'alt'):
            motifs_df[key] = np.where(
                motifs_df['strand'] == '-', 
                complement(motifs_df[key]),
                motifs_df[key]
            )
    
        # Add imbalance data
        self.data_df = variants_df[['logit_es_combined', 'group_id', 'min_fdr']].join(motifs_df)

        # Compute preferred allele
        self.data_df["prefered_allele"] = np.where(
            self.data_df['logit_es_combined'] >= 0,
            self.data_df["ref"],
            self.data_df["alt"])
        self.data_df['ddg'] = self.data_df.eval('ref_score - alt_score')
        
        self.result_df = None

    @staticmethod
    def logtr(arr):
        return 1 - 1 / (1 + np.power(2, np.array(arr)))

    @staticmethod
    def get_concordant(x, y, expected_es=0, x_mar=0):
        diff = x - expected_es
        unmasked_values = (np.abs(diff) >= x_mar) & (y != 0)
        return ((y * diff > 0) & unmasked_values).sum() / unmasked_values.sum()

    @staticmethod
    def set_index(df):
        if len(df.index) == 0:
            return df
        df['variant_id'] = df[['#chr', 'start', 'end', 'ref', 'alt']].astype(str).agg('@'.join, axis=1)
        return df.set_index('variant_id')


    def get_annotations(self, snvs):

        snvs_in = snvs[(snvs['within'] == 1)]

        try:
            if len(snvs_in) == 0:
                raise NoDataException()

            predictor_array = snvs_in['logit_es_combined'].to_numpy()
            expected_es = 0

            if (predictor_array == 0).all() == 0:
                raise NoDataException()
            
            X = predictor_array
            Y = snvs_in['ddg']

            X = sm.add_constant(X)

            model = sm.OLS(Y, X).fit()

            outliers = model.get_influence().cooks_distance[0] > (4/X.shape[0])

            model = sm.OLS(Y[~outliers], X[~outliers]).fit()
            
            r2 = model.rsquared

            concordance = self.get_concordant(
                predictor_array, 
                snvs_in["ddg"], 
                expected_es=expected_es
            )
            
            ref_bias = (predictor_array > 0).sum() / (predictor_array != 0).sum()
        except NoDataException:
            r2, concordance, ref_bias = np.nan, np.nan, np.nan

        return r2, concordance, ref_bias


    #Enrichment code
    def calc_enrichment(self, group_df, imbalanced):
        n_shuffles = 1000
        bins = np.arange(group_df['offset'].min(), group_df['offset'].max() + 2) # add 2 to length: one for '0' and one for last element

        n_all = np.histogram(group_df['offset'], bins=bins)[0]
        n_imbalanced = np.histogram(group_df['offset'][imbalanced], bins=bins)[0]
        
        all_inside = n_all[self.flank_width:-self.flank_width]
        imbalanced_inside = n_imbalanced[self.flank_width:-self.flank_width]

        total_inside = np.nansum(all_inside)
        total_imbalanced_inside = np.nansum(imbalanced_inside)

        log_odds = np.log2(total_imbalanced_inside) - np.log2(total_inside - total_imbalanced_inside) - \
            np.log2(np.nansum(n_imbalanced) - total_imbalanced_inside) + \
                np.log2(np.nansum(n_all) - np.nansum(n_imbalanced) - total_inside + total_imbalanced_inside)
    

        perm = np.zeros(n_shuffles)
        perm_per_nt = np.zeros((n_shuffles, len(bins)-1))

        for i in range(n_shuffles):
            n_exp_imbalanced = np.histogram(group_df['offset'][np.random.permutation(imbalanced)], bins=bins)[0] + 1
            n_exp_not_imbalanced = n_all - n_exp_imbalanced +1 

            perm[i] = np.log2( (n_exp_imbalanced[self.flank_width:-self.flank_width].sum() / n_exp_imbalanced.sum()) / (n_exp_not_imbalanced[self.flank_width:-self.flank_width].sum() / n_exp_not_imbalanced.sum()) )
            perm_per_nt[i,:] = np.log2( (n_exp_imbalanced / np.sum(n_exp_imbalanced)) / (n_exp_not_imbalanced / np.sum(n_exp_not_imbalanced)))

        pval = -stats.norm.logsf(log_odds, loc=np.nanmean(perm, axis=0), scale=np.nanstd(perm, axis=0))/np.log(10)

        # old
        # bins = np.arange(group_df['offset'].min(), group_df['offset'].max() + 2) # add 2 to length: one for '0' and one for last element
        # n_all = np.histogram(group_df['offset'], bins=bins)[0]
        # n_imbalanced = np.histogram(group_df['offset'][imbalanced], bins=bins)[0]

        # all_inside = n_all[self.flank_width:-self.flank_width]
        # imbalanced_inside = n_imbalanced[self.flank_width:-self.flank_width]

        # total_inside = np.nansum(all_inside)
        # total_imbalanced_inside = np.nansum(imbalanced_inside)
        # if total_imbalanced_inside == 0 or total_inside - total_imbalanced_inside == 0:
        #     raise NoDataException()

        # log_odds = np.log2(total_imbalanced_inside) - np.log2(total_inside - total_imbalanced_inside) - \
        #     np.log2(np.nansum(n_imbalanced) - total_imbalanced_inside) + \
        #         np.log2(np.nansum(n_all) - np.nansum(n_imbalanced) - total_inside + total_imbalanced_inside)
        # pval = -stats.hypergeom.logsf(
        #     total_imbalanced_inside,
        #     np.nansum(n_all),
        #     np.nansum(n_imbalanced),
        #     total_inside
        # ) / np.log(10)

        
        return [
            log_odds,
            pval,
            total_inside,
            np.nansum(imbalanced_inside),
            np.nanmedian(all_inside),
            np.nansum(imbalanced_inside >= 7)
        ]


    def get_stats(self, group_df):
        imbalanced_index = group_df['min_fdr'] <= self.fdr_tr
        try:
            if imbalanced_index.sum() == 0:
                raise NoDataException()

            lin_model_stats = self.get_annotations(group_df[imbalanced_index])

            # log_odds, pval, n_inside, n_imb_inside, n_median_inside, n_inside_more_7
            enrichment_stats = self.calc_enrichment(group_df, imbalanced_index)
            data = [
                *enrichment_stats,
                *lin_model_stats
            ]
        except NoDataException:
            data = [np.nan] * (len(self.result_columns) - 2)

        return pd.Series([group_df.name[0], group_df.name[1], *data], self.result_columns)

    def get_motif_stats(self):
        print('Grouping by and applying')
        return self.data_df.groupby(['motif', 'group_id']).progress_apply(self.get_stats)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Converts matrix in txt format to npy')
    parser.add_argument('variants', help='Path to file with CAVs significance')
    parser.add_argument('motifs', help='Path to file with motif hits at variants')
    parser.add_argument('outpath', help='Path to output file')
    parser.add_argument('--flank_width', help='Width of flanks', type=int, default=20)
    parser.add_argument('--fdr', help='FDR threshold for CAVs', type=float, default=0.05)
    args = parser.parse_args()

    data_holder = MotifEnrichment(
        args.variants, args.motifs,
        fdr_tr=args.fdr,
        flank_width=args.flank_width)
    res_df = data_holder.get_motif_stats()
    res_df.to_csv(args.outpath, index=False, sep='\t')