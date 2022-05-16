import argparse
import pandas as pd
from helpers import alleles

def main(input_file, output_file, fdr_tr):
    df = pd.read_table(input_file)
    df = df[df[[f'fdrp_bh_{allele}' for allele in alleles]].min(axis=1) >= fdr_tr]
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate pvalue for model')
    parser.add_argument('-I', help='BABACHI annotated BED file with SNPs')
    parser.add_argument('-O', help='File to save calculated p-value into')
    parser.add_argument('--fdr', type=float, help='Max cover threshold for fdr', default=0.05)
    args = parser.parse_args()
    main(args.I, args.O, args.fdr)