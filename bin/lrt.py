import pandas as pd
import argparse
import numpy as np
from scipy.stats import f
from scipy.special import expit

from aggregation import starting_columns, aggregate_effect_size

# Likelihoods
# L0 <- 'es = 0' model
# L1 <- 'es = mean' model
# L2 <- 'es = mean|group' model

result_columns = [
    *starting_columns,
    'group_id', 'n_groups',
    'log_p_overall',
    'log_p_differential',
    'overall_es', 'group_es', 'group_es_std', 'N',
    'Q0', 'Qtotal', 'Qg',
    'F1', 'F2'
]

class ANOVA:
    def __init__(self, melt, min_indivs_per_group=3, min_groups_per_variant=2, coverage_tr=10):
        self.min_indivs_per_group = min_indivs_per_group
        self.min_groups_per_variant = min_groups_per_variant
        self.coverage_tr = coverage_tr

        melt['variant_id'] = melt['#chr'] + "@" + melt['end'].astype(str) + "@" + melt['alt']       

        testable_pairs = self.find_testable_pairs(melt)
        # filter only testable variants + cell_types
        self.tested_melt = melt.merge(testable_pairs)
        print(f'Testing {self.tested_melt["variant_id"].nunique()} variants')
        if self.tested_melt["variant_id"].nunique() == 0:
            print('No variants for LRT')
        self.tested_melt.drop(columns='variant_id', inplace=True)

    
    def get_testable_snps(self):
        return self.tested_melt

    def test_group(self, df):
        es2 = aggregate_effect_size(df['es'], df['inverse_mse'])
        Qg = np.average(np.square(df['es'] - es2), weights=df['inverse_mse'])
        return pd.Series(
            [es2, np.sqrt(Qg), Qg * df['inverse_mse'].sum()],
            ['group_es', 'group_es_std', 'per_group_Qg']
        )

    def test_snp(self, df):
        n_groups = df['group_id'].nunique()
        N = len(df.index)
        es_mean = aggregate_effect_size(df['es'], df['inverse_mse'])
        Q0 =  np.average(np.square(df['es'] - 0.5), weights=df['inverse_mse'])
        Qtotal =  np.average(np.square(df['es'] - es_mean), weights=df['inverse_mse'])
        Wsum = np.sum(df['inverse_mse'])
        return pd.Series([es_mean, Qtotal, Q0, n_groups, N, Wsum], ['overall_es', 'Qtotal', 'Q0', 'n_groups', 'N', 'Wsum'])

    def find_testable_pairs(self, df):
        # number of samples for particular group with variant_id present
        samples_num = df[['variant_id', 'group_id', 'indiv_id', 'coverage']].groupby(
            ['group_id', 'variant_id']
        ).agg(
            max_coverage=('coverage', 'max'),
            indiv_count=('indiv_id', 'nunique')
        )
        # for each group, variant is present in >= 3 indivs and max_coverage >= 10
        samples_num = samples_num[samples_num.eval(f"max_coverage >= {self.coverage_tr} & indiv_count >= {self.min_indivs}")]
        # of variants for particular group
        groups_per_variant = samples_num.reset_index().value_counts('variant_id')

        # variants present in >=2 groups
        tested_ids = groups_per_variant[groups_per_variant >= self.min_groups_per_variant].index

        testable_variant_group_pairs = samples_num.reset_index()

        return testable_variant_group_pairs[
            testable_variant_group_pairs['variant_id'].isin(tested_ids)
        ]


    def run_anova(self):
        # Total aggregation (find constitutive CAVs)
        res = self.tested_melt.groupby(
            [*starting_columns, 'group_id']
        ).progress_apply(
            self.test_group
        ).reset_index()

        result = self.tested_melt.groupby(starting_columns).progress_apply(
            self.test_snp
        ).join(
            res.groupby(starting_columns).agg(Qg=('per_group_Qg', 'sum'))
        ).reset_index().merge(res)

        result['Qg'] = result['Qg'] / result['Wsum']
    
        assert len(result) == len(res)

        result['F1'] = result.eval('(Q0 - Qtotal) / (Qtotal / (N - 2))')
        result['F2'] = result.eval('(Qtotal - Qg) / (n_groups - 1) / (Qg / (N - n_groups))')
        result['log_p_overall'] = f.logsf(result['F1'], dfn=1, dfd=result['N'] - 2)
        result['log_p_differential'] = f.logsf(
            result['F2'],
            dfn=result['n_groups'] - 1,
            dfd=result['N'] - result['n_groups']
        )
        return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate ANOVA for tested CAVs')
    parser.add_argument('input_data', help='Non-aggregated file with tested CAVs')
    parser.add_argument('metadata', help='Samples metadata with ag_id to indiv_id correspondence')
    parser.add_argument('prefix', help='Prefix to files to save output files into')
    parser.add_argument('--min_indivs_per_group', type=int, help='Number of indivs in each group for the variant', default=3)
    parser.add_argument('--min_groups', type=int, help='Number of groups for the variant', default=2)
    parser.add_argument('--allele_tr', type=int, help='Allelic reads threshold', default=5) # FIXME not used
    parser.add_argument('--chrom', help='Chromosome for parallel execution', default=None)
    parser.add_argument('--coverage_tr', type=int, help='Coverage threshold for at least one variant in the group', default=10)

    args = parser.parse_args()

    input_df = pd.read_table(args.input_data).query('BAD <= 1')
    if args.chrom is not None:
        input_df = input_df[input_df['#chr'] == args.chrom].copy()
    
    ag_id2indiv_id = pd.read_table(args.metadata).set_index('ag_id')['indiv_id'].to_dict()
    input_df['indiv_id'] = input_df['sample_id'].map(ag_id2indiv_id)

    print("Finished reading non-aggregated file, shape:", input_df.shape)
    print("Unique groups:", input_df['group_id'].unique())
    if len(input_df['group_id'].unique()) == 1:
        raise AssertionError("Only one group, LRT is not applicable")

    data_wrapper = ANOVA(
        input_df,
        min_indivs_per_group=args.min_indivs_per_group,
        min_groups_per_variant=args.min_groups,
        coverage_tr=args.coverage_tr
    )
    if data_wrapper.get_testable_snps().empty:
        result = pd.DataFrame([], columns=result_columns)
    else:
        result = data_wrapper.run_anova()[result_columns]
    
    data_wrapper.get_testable_snps().to_csv(
        f"{args.prefix}.tested.bed",
        sep='\t',
        index=False
    )

    result.to_csv(f"{args.prefix}.pvals.bed", sep='\t', index=False)
