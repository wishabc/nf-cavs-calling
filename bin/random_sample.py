import pandas as pd
import argparse
import statsmodels.formula.api as smf
import scipy.stats as st
import numpy as np
from statsmodels.stats.multitest import multipletests


maf_bins_fr = [0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2]
mutations = ['C>A', 'C>T', 'C>G', 'T>A']


def get_sampling_df(df):
    variants_df = df.reset_index()
    variants_df['sample'] = variants_df.groupby('variant_id').cumcount()
    variants_df = variants_df.set_index(['variant_id', 'sample']).sort_index()

    non_unique_variants_ids = variants_df.loc[(slice(None), 1), :].index.get_level_values('variant_id')
    non_unique_n_aggregated = df['variant_id'].value_counts().loc[non_unique_variants_ids]
    non_unique_df_index = variants_df.loc[(non_unique_variants_ids, slice(None)), :].index
    return variants_df, non_unique_n_aggregated, variants_df.index.difference(non_unique_df_index)


def wilson(p, n, z = 1.96):
    denominator = 1 + z ** 2/n
    centre_adjusted_probability = p + z * z / (2 * n)
    adjusted_standard_deviation = np.sqrt((p * (1 - p) + z * z / (4 * n)) / n)
    
    lower_bound = (centre_adjusted_probability - z * adjusted_standard_deviation) / denominator
    upper_bound = (centre_adjusted_probability + z * adjusted_standard_deviation) / denominator
    return (lower_bound, upper_bound)


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


def get_fraction_trend_from_df(df, maf_bins=maf_bins_fr, estimate=False):
    signature = 'signature1'
    #  ctx_label = 'CTX3'
    
    cts = df.pivot_table(index=signature, columns='imbalanced', values='ID', aggfunc='count')
    ctsx = (pd.isna(cts) | (cts <= 2)).iloc[:, 0] + (pd.isna(cts) | (cts <= 2)).iloc[:, 1]
    singular_signatures = ctsx[ctsx].index

    all_data = df[~df[signature].isin(singular_signatures)
                  & df['MAF'].notna() & (df['MAF'] > min(maf_bins))
                  & (df['MAF'] <= max(maf_bins))]
    
    all_data['MAF_rank'] = all_data['MAF'].rank()
    all_data['imbalanced_numeric'] = all_data['imbalanced'].astype(int)
    all_data['maf_bin'] = pd.cut(all_data['MAF'], bins=maf_bins)
    
    frac_reg = all_data.groupby('maf_bin').apply(frac_bins_data).droplevel(1)
    
    if estimate:
        reduced = smf.logit(f"imbalanced_numeric ~ {signature}*pref_orient + FMR + BAD + mut_rates_roulette", data=all_data).fit(maxiter=50)
        full = smf.logit(f"imbalanced_numeric ~ {signature}*pref_orient + FMR + BAD + mut_rates_roulette + MAF_rank", data=all_data).fit(maxiter=50)
        ctxp = smf.logit(f"imbalanced_numeric ~ {signature}*pref_orient", data=all_data).fit(maxiter=50)

        llr = full.llf - reduced.llf
        pval = st.chi2.sf(llr, df=1)
    
        all_data['reduced'] = reduced.predict(all_data)
        all_data['full'] = full.predict(all_data)
        all_data['ctx'] = ctxp.predict(all_data)
    
        gb = all_data.groupby('maf_bin')[['reduced', 'full', 'ctx']].mean()

        frac_reg[['reduced', 'full', 'ctx']] = gb[['reduced', 'full', 'ctx']]
    else:
        llr = 0
        pval = 1
    frac_reg['llr'] = llr
    frac_reg['pval'] = pval

    return frac_reg


def sample_index(n_aggregated, random_state=42):
    rng = np.random.default_rng(seed=random_state)
    sample_ind = pd.Series(
        rng.integers(n_aggregated),
        index=n_aggregated.index,
        name='sample')
    sample_ind = sample_ind.reset_index().rename(columns={'index': 'variant_id'})
    return pd.MultiIndex.from_frame(sample_ind)

def comp(x):
    return {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }.get(x)

def revcomp(s):
    return ''.join([comp(x) for x in s[::-1]])

def palindromic(ref, alt):
    return {ref, alt} == {'A', 'T'} or {ref, alt} == {'G', 'C'}

def get_mutation_stats(row):
    ref = row['ref']
    alt = row['alt']
    assert row['sequence'][20] == ref
    preceding = row['sequence'][:20]
    following = row['sequence'][21:]
    flank_len = len(preceding)
    if palindromic(ref, alt):
        initial_fwd = f'{ref}>{alt}' in mutations
        ref_orient = True
        fwd = initial_fwd  # if cycle dosen't break
        palindromic_res = [True] + [False] * flank_len
        for i in range(flank_len):
            left = preceding[-i - 1]
            right = following[i]
            if not palindromic(left, right):
                fwd = 'T' not in {left, right} and not (left == right == 'G')
                ref_orient = not (fwd ^ initial_fwd)
                break
            palindromic_res[i + 1] = True
        if initial_fwd:
            sub = f'{ref}>{alt}'
        else:
            sub = f'{alt}>{ref}'
        
    else:        
        palindromic_res = [False] * (flank_len + 1)
        for sub, ref_orient, fwd in [
            (f'{ref}>{alt}', True, True),
            (f'{comp(ref)}>{comp(alt)}', True, False),
            (f'{alt}>{ref}', False, True),
            (f'{comp(alt)}>{comp(ref)}', False, False),
        ]:
            if sub in mutations:
                break
            
    if not fwd:
        preceding, following = revcomp(following), revcomp(preceding)
        
    return pd.Series(list(preceding)[-3:] + list(following)[:3] + [sub, fwd, ref_orient] + palindromic_res[:4])


def make_full_df(input_df, context_df):
    input_df['min_pval'] = input_df[['pval_ref', 'pval_alt']].min(axis=1)
    input_df['variant_id'] = input_df['#chr'] + '_' + input_df['end'].astype(str) + '_' + input_df['alt']
    input_df = input_df.merge(context_df)

    input_df['signature1'] = input_df.apply(lambda row: f"{row['-1']}[{row['sub']}]{row['1']}", axis=1)
    input_df['signature2'] = input_df.apply(lambda row: f"{row['-2']}{row['-1']}[{row['sub']}]{row['1']}{row['2']}", axis=1)
    input_df['signature3'] = input_df.apply(lambda row: f"{row['-3']}{row['-2']}{row['-1']}[{row['sub']}]{row['1']}{row['2']}{row['3']}", axis=1)

    input_df[['RAF', 'AAF']] = input_df[['RAF', 'AAF']].replace('.', np.nan).astype(np.float_)
    input_df = input_df[(input_df['MAF'] != '.') & pd.notna(input_df['MAF'])]

    input_df['ref_is_major'] = input_df['RAF'] >= input_df['AAF']
 
    input_df = input_df[input_df.eval('(ref_counts + alt_counts >= 20) & (FMR <= 0.2) & (MAF != 0)')]
    input_df['es_maj'] = np.where(input_df['ref_is_major'], input_df['es'], -input_df['es'])
    input_df['imbalanced'] = (input_df['FDR'] < 0.05) & (np.abs(input_df['es']) >= np.log2(1.5))
    input_df['MAF_quantile'] = pd.cut(input_df['MAF'], input_df['MAF'].quantile(np.linspace(0, 1, 10)), include_lowest=True)
    input_df['cover'] = input_df.eval('ref_counts + alt_counts')
    input_df['negative'] = (input_df['cover'] >= 30) & (input_df['min_pval'] >= 0.3)
    input_df['pval_quantile'] = pd.cut(input_df['min_pval'], input_df['min_pval'].quantile(np.linspace(0, 1, 10)), include_lowest=True)
    input_df['imbalanced_side_is_major'] = input_df['es_maj'] > 0
    input_df.loc[input_df['es_maj'] == 0, 'imbalanced_side_is_major'] = np.nan
    input_df['revMAF'] = 0.5 - input_df['MAF']

    input_df[
    ['-3', '-2', '-1', '1', '2', '3', 'sub', 'fwd', 'ref_orient', 'palindromic']
        + [f'palindromic_{i}' for i in range(1, 4)]] = input_df.progress_apply(
        get_mutation_stats, axis=1
    )

    return input_df


def main(nonaggregated_df, seed_start=20, seed_step=10):
    # sampling_df - df with 2-level index: [variant_id, count (0-based)]
    print('Making sampling df')
    sampling_df, non_unique_n_aggregated, unique_index = get_sampling_df(nonaggregated_df)

    frac_regs = []
    for seed in range(seed_start, seed_start + seed_step + 1):
        print(f'Processing seed: {seed}')
        sampled_variants_index = sample_index(
            non_unique_n_aggregated,
            seed
        )
        sample_df = sampling_df.loc[sampled_variants_index.union(unique_index)]
        sample_df['FDR'] = multipletests(
                    sample_df['min_pval'],
                    method='fdr_bh'
                )[1]
        sample_df['imbalanced'] = (sample_df['FDR'] <= 0.05) #& (np.abs(sample_df['es']) >= np.log2(1.5))
        frac = get_fraction_trend_from_df(sample_df, estimate=True)
        frac_regs.append(frac)

    return pd.concat(frac_regs)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Sampling one of recurrent variants")
    parser.add_argument('-I', help='Non-aggregated BED file')
    parser.add_argument('-c', help='File with context annotation')
    parser.add_argument('-O', help='File to save calculated metrics')
    parser.add_argument('--start', type=int, help='Start value for seed', default=10)
    parser.add_argument('--step', type=int, help='Step size for seed values', default=10)

    args = parser.parse_args()
    input_df = pd.read_table(args.I)
    context_df = pd.read_table(args.c, header=None, names=['#chr', 'start', 'end', 'sequence'])
    
    input_df = make_full_df(input_df, context_df)
    df = main(input_df, seed_start=args.start, seed_step=args.step)
    df.to_csv(args.O, sep='\t', index=False)