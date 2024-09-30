import pandas as pd
import sys


def main(bed_file):
    return bed_file[bed_file.eval('alt == alt_mr')]

if __name__ == '__main__':
    bed_file = pd.read_table(sys.argv[1])
    filt_bed = main(bed_file)
    filt_bed[['#chr', 'start', 'end', 'ID', 'ref',
        'alt', 'mut_rates_roulette', 'mut_rates_gnomad']
        ].to_csv(sys.argv[2], sep='\t', index=False)
