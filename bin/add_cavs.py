import pandas as pd
import argparse
from helpers import starting_columns


def set_index(df):
    if not df.empty:
        df.index = df.apply(lambda row: "@".join(map(str, [row[x] for x in starting_columns + ['sample_id']])), axis=1)
    return df

def main(new_badmap, old_badmap, output):
    old_df = set_index(pd.read_table(old_badmap))
    if new_badmap is None:
        old_badmap.to_csv(output, sep='\t', index=False)
        return
    new_df = set_index(pd.read_table(new_badmap))
    if new_df.empty:
        new_df.to_csv(output, sep='\t', index=False)
        return
    imputed_cavs = old_df.loc[old_df.index.difference(new_df.index)]
    df = pd.concat([new_df, imputed_cavs])
    if len(df.index) != len(old_df.index):
        print(len(df.index), len(old_df.index), len(old_df.index.difference(new_df.index)),
         len(new_df.index), len(imputed_cavs.index))
        raise AssertionError
    df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge results of two BADmaps calling iterations')
    parser.add_argument('-n', help='New BAD intersection file', default=None)
    parser.add_argument('-o', help='Old BAD intersection file')
    parser.add_argument('--output', help='Output file name')
    args = parser.parse_args()
    main(new_badmap=args.n, old_badmap=args.o, output=args.output)