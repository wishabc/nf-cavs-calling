import pandas as pd
import argparse
from estimate_mse import calc_mse
import numpy as np

def read_indicator(file_path):
    return np.loadtxt(file_path, dtype=str)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add annotations to the bed file.')
    parser.add_argument('bed_file', type=str, help='Path to the bed file')
    parser.add_argument('output', type=str, help='Path to the output file')
    parser.add_argument('--hotspots', type=str, help='Path to the in-hotspot indicator file')
    parser.add_argument('--footprints', type=str, help='Path to the in-footprint indicator file')
    args = parser.parse_args()

    print('Reading files')

    df = pd.read_table(args.bed_file)
    initial_df_len = len(df.index)
    hotspots = read_indicator(args.hotspots)
    footprints = read_indicator(args.footprints)
    assert len(hotspots) == len(footprints) == initial_df_len, f"Hotspots {len(hotspots)} and footprints {len(footprints)} should have the same length as the bed file {initial_df_len}."
    df['hotspots'] = hotspots
    df['footprints'] = footprints
    
    data = df[['BAD', 'coverage']].drop_duplicates()
    data['inverse_mse'] = data.progress_apply(lambda row: 1 / calc_mse(row['coverage'], row['BAD']).mean(), axis=1)
    df = df.drop(columns=['inverse_mse', 'inverse_mse_x', 'inverse_mse_y'], errors='ignore').merge(data, on=['BAD', 'coverage'])
    assert initial_df_len == len(df.index)
    df.to_csv(args.output, index=False, sep='\t')
    