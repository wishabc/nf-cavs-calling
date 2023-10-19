import pandas as pd
import argparse


def main(df, fdr_tr=0.1):
    tested_vars = df['group_id'].value_counts().rename('tested')
    imbalanced_vars = df[df['min_fdr'] <= fdr_tr]['group_id'].value_counts().rename('imbalanced')
    result = pd.concat([tested_vars, imbalanced_vars], axis=1).reset_index(names='group_id')
    result['percentage'] = result.eval('imbalanced / tested')
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Collect CAVs statistics')
    parser.add_argument('I', help='File with aggregated pvalues')
    parser.add_argument('O', help='File to save statistics into')
    parser.add_argument('--fdr',
        help='FDR threshold for variant to be considered a CAV',
        type=float, default=0.1)
    args = parser.parse_args()
    aggregated_variants = pd.read_table(args.I)
    res_df = main(aggregated_variants[aggregated_variants['min_fdr'].notna()], args.fdr)
    res_df.to_csv(args.O, sep='\t', index=False)
