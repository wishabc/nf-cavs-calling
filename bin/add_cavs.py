import pandas as pd
import argparse
from aggregation import starting_columns


def main(new_badmap, old_badmap, output):
    group_cols = [*starting_columns, 'sample_id']
    old_df = pd.read_table(old_badmap, low_memory=False).set_index(group_cols)
    print(old_df.index)
    if new_badmap is None:
        old_df.to_csv(output, sep='\t', index=False)
        return
    new_df = pd.read_table(new_badmap, low_memory=False).set_index(group_cols)
    if new_df.empty:
        old_df.to_csv(output, sep='\t', index=False)
        return
    assert len(new_df.index.difference(old_df.index)) == 0
    old_df.loc[new_df.index] = new_df
    old_df.reset_index().to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Merge results of two BADmaps calling iterations')
    parser.add_argument('-n', help='New BAD intersection file', default=None)
    parser.add_argument('-o', help='Old BAD intersection file')
    parser.add_argument('--output', help='Output file name')
    args = parser.parse_args()
    main(new_badmap=args.n, old_badmap=args.o, output=args.output)