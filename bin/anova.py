import os
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import pandas as pd
import argparse


def read_non_aggregated_files(dir_path):
    tables = []
    for file in os.listdir(dir_path):
        group_id = file.split('.')[0]
        tb_df = pd.read_table(os.path.join(dir_path, file))
        tb_df['group_id'] = group_id
        tables.append(tb_df)
    res = pd.concat(tables)
    res['variant_id'] = res['#chr'] + '_' + res['end'].astype(str) + '_' + res['alt']
    return res

def linear_reg(df):
    snp_data = df.iloc[0]
    return pd.DataFrame({
        **{x: [snp_data[x]] for x in ['#chr', 'start', 'end', 'ID', 'ref', 'alt']},
        'groups': '|'.join(df['group_id'].tolist()),
        'anova_pvalue': [smf.ols('es ~ C(group_id)', data=df).fit().f_pvalue],
    })


def find_testable_pairs(df, min_samples, min_groups_per_variant):
    # # of samples for particular group with variant_id present
    samples_num = df.value_counts(['group_id', 'variant_id'])
    # for each group, variant is present in >= 3 samples
    samples_num = samples_num[samples_num >= min_samples]
    
    # of variants for particular group
    groups_per_variant = samples_num.reset_index().value_counts('variant_id')

    # variants present in >=2 groups
    tested_ids = groups_per_variant[groups_per_variant >= min_groups_per_variant].index

    testable_variant_group_pairs = samples_num.reset_index()

    return testable_variant_group_pairs[
         testable_variant_group_pairs['variant_id'].isin(tested_ids)
    ]


def main(melt, min_samples=3, min_groups=2, cover_tr=20):
    melt = melt[melt.eval(f'ref_counts + alt_counts >= {cover_tr}')]
    testable_pairs = find_testable_pairs(melt, min_samples=min_samples,
        min_groups_per_variant=min_groups)
    # filter only testable variants + cell_types
    tested_melt = melt.merge(
        testable_pairs, on=['variant_id', 'group_id']
    )
    return tested_melt.groupby('variant_id').apply(linear_reg)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate ANOVA for tested CAVs')
    parser.add_argument('input_file', help='Non-aggregated file with tested CAVs')
    parser.add_argument('output_file', help='File to save calculated ANOVA p-value into')
    parser.add_argument('--ct', type=int, help='Cover threshold for fdr', default=20)
    parser.add_argument('--min_samples', type=int, help='Number of samples in each group for the variant', default=3)
    parser.add_argument('--min_groups', type=int, help='Number of groups for the variant', default=2)
    args = parser.parse_args()

    melt = read_non_aggregated_files(args.input_file)
    min_samples = args.min_samples # of samples
    min_groups_per_variant = args.min_groups
    fdr_cov_tr = args.ct
    anova_results = main(melt, cover_tr=fdr_cov_tr,
        min_samples=min_samples, min_groups=min_groups_per_variant)
    
    anova_results['anova_fdr'] = multipletests(anova_results['anova_pvalue'],
        method='fdr_bh')[1]
    anova_results.to_csv(args.output_file, sep='\t', index=False)
