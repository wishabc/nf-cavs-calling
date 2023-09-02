import pandas as pd
import argparse
import statsmodels.formula.api as smf
import scipy.stats as st
import numpy as np
from tqdm import tqdm
from aggregation import calc_fdr, starting_columns

tqdm.pandas()

maf_bins_fr = [0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2]


def get_sampling_df(df):
    variants_df = df.reset_index().copy()
    variants_df['sample'] = variants_df.groupby(starting_columns).cumcount()
    variants_df = variants_df.set_index(['variant_id', 'sample']).sort_index()

    non_unique_variants_ids = variants_df.loc[(slice(None), 1), :].index.get_level_values('variant_id')
    non_unique_n_aggregated = df['variant_id'].value_counts().loc[non_unique_variants_ids]
    non_unique_df_index = variants_df.loc[(non_unique_variants_ids, slice(None)), :].index
    return variants_df, non_unique_n_aggregated, variants_df.index.difference(non_unique_df_index)


def wilson(p, n, z=1.96):
    denominator = 1 + z ** 2/n
    centre_adjusted_probability = p + z * z / (2 * n)
    adjusted_standard_deviation = np.sqrt((p * (1 - p) + z * z / (4 * n)) / n)
    
    lower_bound = (centre_adjusted_probability - z * adjusted_standard_deviation) / denominator
    upper_bound = (centre_adjusted_probability + z * adjusted_standard_deviation) / denominator
    return lower_bound, upper_bound


def frac_bins_data(df):
    pos = (df['imbalanced']).sum()
    n = len(df.index)
    p = pos / n
    if n == 0:
        lb, ub = 0, 1
    else:
        lb, ub = wilson(p, n)
    return pd.DataFrame({
        'imb_frac': [p],
        'frac_lower_bound': [lb],
        'frac_upper_bound': [ub],
        'num_observations': [n],
    })


def get_fraction_trend_from_df(all_data):

    # reduced = smf.logit(f"imbalanced_numeric ~ signature1 * pref_orient + FMR + BAD + mut_rates_roulette", data=all_data).fit(maxiter=50)
    # full = smf.logit(f"imbalanced_numeric ~ signature1 * pref_orient + FMR + BAD + mut_rates_roulette + MAF_rank", data=all_data).fit(maxiter=50)
    # ctxp = smf.logit(f"imbalanced_numeric ~ signature1 * pref_orient", data=all_data).fit(maxiter=50)
    
    # all_data['reduced'] = reduced.predict(all_data)
    # all_data['full'] = full.predict(all_data)
    # all_data['ctx'] = ctxp.predict(all_data)
    
    # gb = all_data.groupby('maf_bin')[['reduced', 'full', 'ctx']].mean()

    frac_reg = all_data.groupby('maf_bin').apply(frac_bins_data).reset_index()
    # frac_reg[['reduced', 'full', 'ctx']] = gb[['reduced', 'full', 'ctx']]
    # frac_reg['llr'] = full.llf - reduced.llf
    # frac_reg['pval'] = st.chi2.sf(frac_reg['llr'], df=1)

    return frac_reg


def sample_index(n_aggregated, random_state=42):
    rng = np.random.default_rng(seed=random_state)
    sample_ind = pd.Series(
        rng.integers(n_aggregated),
        index=n_aggregated.index,
        name='sample')
    sample_ind = sample_ind.reset_index().rename(columns={'index': 'variant_id'})
    return pd.MultiIndex.from_frame(sample_ind)


def main(input_df, annotation_df, seed_start=20, seed_step=10):
    # sampling_df - df with 2-level index: [variant_id, count (0-based)]
    # This script can be further optimized by moving preprocessing
    # to a separate script
    print('Preprocessing df')
    input_df = input_df[input_df['is_tested']].copy()

    input_df[['RAF', 'AAF']] = input_df[['RAF', 'AAF']].apply(
        lambda x: pd.to_numeric(x, errors='coerce')
    )
    input_df.dropna(subset=['RAF', 'AAF'], inplace=True)
    input_df['MAF'] = input_df[['RAF', 'AAF']].min(axis=1, skina=False)
    input_df['MAF_rank'] = input_df['MAF'].rank()
    input_df['maf_bin'] = pd.cut(input_df['MAF'], bins=maf_bins_fr)

    es_mean = input_df.groupby(starting_columns)['es'].mean().reset_index().rename(
        columns={'es': 'es_weighted_mean'}
    )
    input_df = input_df.merge(annotation_df).merge(es_mean)
    input_df['pref_orient'] = np.where(
        input_df['ref_orient'], 
        input_df['es_weighted_mean'] > 0,
        input_df['es_weighted_mean'] < 0
    )
    input_df['variant_id'] = input_df[starting_columns].astype(str).agg('@'.join, axis=1)
    
    cts = input_df.pivot_table(
        index='signature1', 
        columns='imbalanced',
        values='ID',
        aggfunc='count'
    )
    
    ctsx = (pd.isna(cts) | (cts <= 2)).iloc[:, 0] + (pd.isna(cts) | (cts <= 2)).iloc[:, 1]
    singular_signatures = ctsx[ctsx].index

    input_df = input_df[~input_df['signature1'].isin(singular_signatures) & ~pd.isna(input_df['maf_bin'])].copy()
    

    sampling_df, non_unique_n_aggregated, unique_index = get_sampling_df(input_df)
    print('Preprocessing finished')

    frac_regs = []

    for seed in tqdm(list(range(seed_start, seed_start + seed_step))):
        print(f'Processing seed: {seed}')
        sampled_variants_index = sample_index(
            non_unique_n_aggregated,
            seed
        )
        sample_df = sampling_df.loc[sampled_variants_index.union(unique_index)]
        sample_df = calc_fdr(sample_df, prefix='pval_')
        sample_df['imbalanced'] = sample_df['min_fdr'] <= 0.05
        sample_df['imbalanced_numeric'] = sample_df['imbalanced'].astype(int)
    
        frac = get_fraction_trend_from_df(sample_df)
        frac_regs.append(frac)

    return pd.concat(frac_regs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Sampling one of recurrent variants")
    parser.add_argument('I', help='Non-aggregated BED file')
    parser.add_argument('a', help='File with annotations')
    parser.add_argument('O', help='File to save calculated metrics')
    parser.add_argument('--noncpg', help='Use only non-cpg', default=False, action="store_true")
    parser.add_argument('--start', type=int, help='Start value for seed', default=10)
    parser.add_argument('--step', type=int, help='Step size for seed values', default=10)

    args = parser.parse_args()

    input_df = pd.read_table(args.I)
    annotation_df = pd.read_table(args.a)
    df = main(input_df, annotation_df, seed_start=args.start, seed_step=args.step)
    df.to_csv(args.O, sep='\t', index=False)