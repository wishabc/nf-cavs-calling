import argparse


def main():
    pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate ANOVA for tested CAVs')
    parser.add_argument('tested_variants', help='Non-aggregated file with tested CAVs')
    parser.add_argument('prefix', help='Prefix to files to save output files into')
    parser.add_argument('--min_samples', type=int, help='Number of samples in each group for the variant', default=3)
    parser.add_argument('--min_groups', type=int, help='Number of groups for the variant', default=2)
    parser.add_argument('--allele_tr', type=int, help='Allelic reads threshold', default=5)
    parser.add_argument('--chrom', help='Chromosome for parallel execution', default=None)

    args = parser.parse_args()