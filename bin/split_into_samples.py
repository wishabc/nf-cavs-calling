import pandas as pd
import sys


def save_as_tsv(df, suffix):
    df.to_csv(f"{df.name}.{suffix}", index=False, sep='\t')

def main(df, suffix):
    df.groupby('sample_id').apply(lambda x: save_as_tsv(x, suffix))

if __name__ == '__main__':
    df = pd.read_table(sys.argv[1])
    suf = sys.argv[2]
    main(df, suf)
   