import pandas as pd
import argparse
import numpy as np
import scipy.stats as st
from scipy.special import expit
from aggregation import calc_fdr_pd, check_if_tested
from AI_statistics.vectorized_estimators import logit_es, es_estimate_vectorized, estimate_w_null, calc_bimodal_pvalues
from tqdm import tqdm

tqdm.pandas()


updated_columns = ['coverage', 'w', 'es', 'logit_es', 'pval_ref', 'pval_alt', 'min_pval', 'FDR_sample']


def main(df: pd.DataFrame, coverage_tr=15):
    df = df.drop(columns=updated_columns, errors='ignore')
    # Check if empty
    result_columns = [*df.columns, *updated_columns]
    if df.empty:
        return pd.DataFrame([], columns=result_columns)

    df['coverage'] = df.eval('ref_counts + alt_counts')
    df = df.groupby('sample_id').progress_apply(check_if_tested, max_cover_tr=coverage_tr).reset_index(drop=True)
    
    df['w'] = estimate_w_null(df['ref_counts'], df['coverage'], df['BAD'])
    df['es'] = es_estimate_vectorized(df['ref_counts'], df['coverage'], df['BAD'], df['w'])
    df['logit_es'] = logit_es(df['es'])

    log_bad = np.log(df['BAD'])

    dist1 = st.binom(df['coverage'], expit(log_bad))
    dist2 = st.binom(df['coverage'], expit(-log_bad))
    log_pval_ref, log_pval_alt, log_pval_both = calc_bimodal_pvalues(
        dist1,
        dist2,
        df['ref_counts'], 
        df['w']
    )
    df = df.assign(
        pval_ref=np.exp(log_pval_ref),
        pval_alt=np.exp(log_pval_alt),
        min_pval=np.exp(log_pval_both),
    )
    df['FDR_sample'] = df.groupby('sample_id')['min_pval'].transform(calc_fdr_pd)
    return df.reset_index(drop=True)[result_columns]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('--coverage_threhold', type=int, help="""Coverage threshold for variants to calculate per-sample q-values.
                        Expected to be a positive integer""", default=15)
    args = parser.parse_args()
    input_df = pd.read_table(args.I, low_memory=False)
    modified_df = main(
        input_df,
        coverage_tr=args.coverage_threhold
    )
    modified_df.to_csv(args.O, sep='\t', index=None)