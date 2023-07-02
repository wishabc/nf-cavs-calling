import pandas as pd
import argparse
from helpers import starting_columns


def set_index(df):
    if not df.empty:
        df.index = df.apply(lambda row: "@".join(map(str, [row[x] for x in [*starting_columns, 'sample_id']])), axis=1)
    return df

def main(new_badmap, old_badmap, output):
    old_df = set_index(pd.read_table(old_badmap, low_memory=False))
    if new_badmap is None:
        old_df.to_csv(output, sep='\t', index=False)
        return
    new_df = set_index(pd.read_table(new_badmap, low_memory=False))
    if new_df.empty:
        new_df.to_csv(output, sep='\t', index=False)
        return
    assert len(new_df.index.difference(old_df.index)) == 0
    old_df.loc[new_df.index] = new_df
    old_df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge results of two BADmaps calling iterations')
    parser.add_argument('-n', help='New BAD intersection file', default=None)
    parser.add_argument('-o', help='Old BAD intersection file')
    parser.add_argument('--output', help='Output file name')
    args = parser.parse_args()
    main(new_badmap=args.n, old_badmap=args.o, output=args.output)