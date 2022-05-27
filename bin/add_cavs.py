import pandas as pd
import argparse

def add_key(df):
    df['key'] = df['#chr'] + '@' + df['start']
    return df

def main(new_badmap, old_badmap, output):
    new_df = add_key(pd.read_table(new_badmap))
    old_df = add_key(pd.read_table(old_badmap))
    rs_ids = new_df['key'].unique()
    
    imputed_cavs = old_df[~old_df['key'].isin(rs_ids)]
    df = pd.concat([new_df, imputed_cavs])
    if len(df.index) != len(old_df.index):
        print(len(df.index), len(old_df.index), len(rs_ids), len(imputed_cavs.index))
        raise AssertionError
    df[[x for x in df.columns if x != 'key']].to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add CAVs')
    parser.add_argument('-n', help='New BAD map file')
    parser.add_argument('-o', help='Old BAD map file')
    parser.add_argument('--output', help='Output file name')
    args = parser.parse_args()
    main(new_badmap=args.n, old_badmap=args.o, output=args.output)