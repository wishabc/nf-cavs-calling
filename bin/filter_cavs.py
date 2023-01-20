import argparse
import pandas as pd
from helpers import starting_columns

def main(agg_file, snps_file, output_file, fdr_tr):
    df = pd.read_table(agg_file)
    if not df.empty:
        df = df[(df['min_fdr'] >= fdr_tr) | pd.isna(df['min_fdr'])]
    snps_df = pd.read_table(snps_file)
    if not snps_df.empty:
        snps_df = snps_df.merge(df[starting_columns], how='inner')[starting_columns + ['ref_counts', 'alt_counts', 'sample_id', 'AAF', 'RAF', 'FMR']]

    snps_df.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-a', help='Aggregated CAVs file')
    parser.add_argument('-b', help='BABACHI annotated SNPs BED file')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('--fdr', type=float, help='FDR threshold for CAVs', default=0.05)
    args = parser.parse_args()
    main(args.a, args.b, args.O, args.fdr)