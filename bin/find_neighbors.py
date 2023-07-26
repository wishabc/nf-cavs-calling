import pandas as pd
import sys

columns = ['chr', 'pair_first_start', 'pair_first', 'pair_second', 'distance', 'start_es', 'end_es', 'sample_id']

# assumes df is sorted
def get_pairs(df, group_name):
    # returns indeces of the pairs in dataframe + distance
    combs = list(zip(range(len(df.index) - 1), range(1, len(df.index))))
    for comb in combs:
        assert comb[0] < comb[1]
    if len(combs) == 0:
        return pd.DataFrame([], columns=columns)
    starts, ends = map(list, zip(*combs))
    positions = df.end.to_numpy()
    es = df.es.to_numpy()
    result = [
        group_name[0],
        (positions - 1)[starts], positions[starts],
        positions[ends], positions[ends] - positions[starts],
        es[starts], es[ends], 
        group_name[1]
    ]
    
    return pd.DataFrame({y: x for x, y in zip(result, columns)})




def main(sample_df):
    if sample_df.empty:
        return pd.DataFrame([], columns=columns)
    groups = sample_df.groupby(['#chr', 'sample_id'])
    pairs = []
    for group_name in list(groups.groups):
        p = get_pairs(groups.get_group(group_name).reset_index(drop=True), group_name)
        if len(p.index) > 0:
            pairs.append(p)
    if len(pairs) == 0:
        return pd.DataFrame([], columns=columns)
    return pd.concat(pairs).reset_index(drop=True)

if __name__ == '__main__':
    sample_df = pd.read_table(sys.argv[1])
    result_df = main(sample_df[sample_df['is_tested']])
    result_df.to_csv(sys.stdout, index=False, sep='\t', header=False)