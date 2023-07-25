import pandas as pd
import argparse
import numpy as np
from scipy.optimize import minimize
from scipy.stats import binom, chi2
from statsmodels.stats.multitest import multipletests
from aggregation import aggregate_pvalues_df, calc_fdr, starting_columns

# Likelihoods
# L0 <- 'es = 0' model
# L1 <- 'es = mean' model
# L2 <- 'es = mean|group' model


result_columns = [
    *starting_columns,
    'group_id',
    'p_overall', 'min_fdr_overall',
    'p_differential', 'differential_FDR',
    'min_fdr_group',
    'es1', 'es2', 'es2_std',
    'DL1', 'DL2'
]
class LRT:
    def __init__(self, melt, allele_tr=5, min_samples=3, min_groups_per_variant=2):
        self.min_samples = min_samples
        self.min_groups_per_variant = min_groups_per_variant
        self.allele_tr = allele_tr

        melt['variant_id'] = melt['#chr'] + '_' + melt['end'].astype(str) + '_' + melt['alt']
        melt['n'] = melt.eval('ref_counts + alt_counts')

        melt['x'] = np.round(
            np.where(
                melt['BAD'] == 1, 
                melt['ref_counts'],
                ((1 - 1 / (1 + np.power(2, melt['es']))) * melt['n'])
            )
        )

        testable_pairs = self.find_testable_pairs(melt)
        # filter only testable variants + cell_types
        self.tested_melt = melt.merge(
            testable_pairs, on=['variant_id', 'group_id']
        )
        print(f'Testing {self.tested_melt["variant_id"].nunique()} variants')
        if self.tested_melt["variant_id"].nunique() == 0:
            print('No variants for LRT')
        
    def get_testable_snps(self):
        return self.tested_melt

    def test_group(self, df):
        alpha = np.log(2)/2
        es2, per_group_L2 = self.get_ml_es_estimation(df['x'], df['n'])
        es2_std = np.cosh(alpha * es2) / (alpha * np.sqrt(df['n'].sum()))
        return pd.Series(
            [df.name, es2, es2_std, per_group_L2],
            ['group_id', 'es2', 'es2_std', 'per_group_L2']
        )

    def test_snp(self, df):
        res = df.groupby('group_id').apply(self.test_group)

        m = df['group_id'].nunique()
        x = df['x'].to_numpy()
        n = df['n'].to_numpy()
        L0 = self.censored_binomial_likelihood(x, n, 0.5).sum()
        L2 = res['per_group_L2'].sum()
        es1, L1 = self.get_ml_es_estimation(x, n)
        res = res.assign(
            variant_id=df.name,
            es1=es1,
            DL1=L1-L0,
            DL2=L2-L1,
            n_groups=m
        )
        return res
    

    def censored_binomial_likelihood(self, xs, ns, p):
        b = binom(ns, p)
        return b.logpmf(xs) - np.log(b.cdf(ns - self.allele_tr) - b.cdf(self.allele_tr - 1))


    def get_ml_es_estimation(self, x, n):
        def target(p):
            return -self.censored_binomial_likelihood(x, n, p).sum()
        p = minimize(target, 0.5, bounds=((0.001, 0.999),)).x[0]
        e = np.log2(p) - np.log2(1 - p)
        return e, -target(p)


    def find_testable_pairs(self, df):
        # # of samples for particular group with variant_id present
        samples_num = df.value_counts(['group_id', 'variant_id'])
        # for each group, variant is present in >= 3 samples
        samples_num = samples_num[samples_num >= self.min_samples]
        
        # of variants for particular group
        groups_per_variant = samples_num.reset_index().value_counts('variant_id')

        # variants present in >=2 groups
        tested_ids = groups_per_variant[groups_per_variant >= self.min_groups_per_variant].index

        testable_variant_group_pairs = samples_num.reset_index().drop(columns=0)

        return testable_variant_group_pairs[
            testable_variant_group_pairs['variant_id'].isin(tested_ids)
        ]


    def run_anova(self):
        # FIXME

        # Total aggregation (find constitutive CAVs)

        constitutive_df = calc_fdr(
            aggregate_pvalues_df(self.tested_melt)
        ).rename(
            columns={'min_fdr': 'min_fdr_overall'}
        )
        # merge with tested variants
        tested_melt = self.tested_melt.merge(
            constitutive_df[[*starting_columns, 'min_fdr_overall']], 
            how='left'
        )

        # LRT (ANOVA-like)
        result = self.tested_melt[['x', 'n', 'variant_id', 'group_id']].groupby(
            'variant_id'
        ).progress_apply(self.test_snp)

        result = tested_melt.merge(result)
        result['p_overall'] = chi2.logsf(result['DL1'], 1)
        result['p_differential'] = chi2.logsf(
            result['DL2'],
            result['n_groups'] - 1
        )

        result['differential_FDR'] = multipletests(
            np.exp(result['p_differential']),
            method='fdr_bh'
        )[1]
        print(len(result.index))

        # set default inividual fdr and find differential snps
        differential_idxs = result['differential_FDR'] <= 0.05
        
        # Group-wise aggregation
        group_wise_aggregation = calc_fdr(
            result[differential_idxs].groupby('group_id').progress_apply(
                aggregate_pvalues_df
            )
        ).rename(
            columns={'min_fdr': 'min_fdr_group'}
        ).reset_index()[[*starting_columns, 'group_id', 'min_fdr_group']]

        return result.merge(group_wise_aggregation, how='left')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate ANOVA for tested CAVs')
    parser.add_argument('input_data', help='Non-aggregated file with tested CAVs')
    parser.add_argument('prefix', help='Prefix to files to save output files into')
    parser.add_argument('--min_samples', type=int, help='Number of samples in each group for the variant', default=3)
    parser.add_argument('--min_groups', type=int, help='Number of groups for the variant', default=2)
    parser.add_argument('--allele_tr', type=int, help='Allelic reads threshold', default=5)
    parser.add_argument('--chrom', help='Chromosome for parallel execution', default=None)

    args = parser.parse_args()

    input_df = pd.read_table(args.input_data)
    input_df = input_df[(input_df['is_tested']) & (True if args.chrom is None else input_df['#chr'] == args.chrom)].copy()
    data_wrapper = LRT(
        input_df,
        min_samples=args.min_samples,
        min_groups_per_variant=args.min_groups,
        allele_tr=args.allele_tr
    )
    if LRT.get_testable_snps().empty:
        result = pd.DataFrame([], columns=result_columns)
    else:
        result = data_wrapper.run_anova()
    
    LRT.get_testable_snps().to_csv(f"{args.prefix}.tested.bed", sep='\t', index=False)
    result.to_csv(f"{args.prefix}.pvals.bed", sep='\t', index=False)
