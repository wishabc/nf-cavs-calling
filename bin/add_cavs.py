import pandas as pd
import argparse

def add_key(df):
    if not df.empty:
        df['key'] = df.apply(lambda row: f"{row['#chr']}@{row['start']}", axis=1)
    return df

def main(new_badmap, old_badmap, output):
    new_df = add_key(pd.read_table(new_badmap))
    if new_df.empty:
        new_df.to_csv(output, sep='\t', index=False)
        return
    old_df = add_key(pd.read_table(old_badmap))
    keys = new_df['key'].unique()
    
    imputed_cavs = old_df[~old_df['key'].isin(keys)]
    df = pd.concat([new_df, imputed_cavs])
    if len(df.index) != len(old_df.index):
        print(len(df.index), len(old_df.index), len(keys), len(imputed_cavs.index))
        raise AssertionError
    df[[x for x in df.columns if x != 'key']].to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add CAVs')
    parser.add_argument('-n', help='New BAD map file')
    parser.add_argument('-o', help='Old BAD map file')
    parser.add_argument('--output', help='Output file name')
    args = parser.parse_args()
    main(new_badmap=args.n, old_badmap=args.o, output=args.output)