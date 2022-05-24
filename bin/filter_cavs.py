import argparse
import pandas as pd
from helpers import alleles

def main(agg_file, snps_file, output_file, fdr_tr):
    df = pd.read_table(agg_file)
    join_columns = ['#chr', 'start', 'end', 'ID', 'ref', 'alt']
    if not df.empty:
        df = df[df[[f'fdrp_bh_{allele}' for allele in alleles]].min(axis=1) >= fdr_tr]
    snps_df = pd.read_table(snps_file)
    if not snps_df.empty:
        res_df = snps_df.join(df, on=join_columns, how='inner')[join_columns + ['ref_counts', 'alt_counts']]
    res_df.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-a', help='Aggregated CAVs file')
    parser.add_argument('-b', help='BABACHI annotated SNPs BED file')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('--fdr', type=float, help='Max cover threshold for fdr', default=0.05)
    args = parser.parse_args()
    main(args.a, args.b, args.O, args.fdr)