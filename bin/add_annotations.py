import pandas as pd
import argparse
from estimate_mse import calc_mse
import numpy as np

def read_indicator(file_path):
    return np.loadtxt(file_path, dtype=int).astype(bool)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add annotations to the bed file.')
    parser.add_argument('bed_file', type=str, help='Path to the bed file')
    parser.add_argument('output', type=str, help='Path to the output file')
    parser.add_argument('--hotspots', type=str, help='Path to the in-hotspot indicator ifile')
    parser.add_argument('--footprints', type=str, help='Path to the in-footprint indicator file')
    args = parser.parse_args()

    df = pd.read_table(args.bed_file)
    initial_df_len = len(df.index)
    hotspots = read_indicator(args.hotspots),
    footprints = read_indicator(args.footprints)
    assert len(hotspots) == len(footprints) == initial_df_len, f"Hotspots {len(hotspots)} and footprints {len(footprints)} should have the same length as the bed file {initial_df_len}."
    df['hotspots'] = hotspots
    df['footprints'] = footprints
    
    data = df[['BAD', 'coverage']].drop_duplicates()
    data['inverse_mse'] = data.progress_apply(lambda row: 1 / calc_mse(row['coverage'], row['BAD']).mean(), axis=1)
    df = df.merge(data, on=['BAD', 'coverage'])
    assert initial_df_len == len(df.index)
    