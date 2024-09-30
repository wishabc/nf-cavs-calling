import pandas as pd
import argparse
from AI_statistics.vectorized_estimators import calc_binom_variance
import numpy as np
from tqdm import tqdm


tqdm.pandas()


def read_indicator(file_path):
    return np.loadtxt(file_path, dtype=str)


def calc_inverse_mse(row):
    n, B = row['coverage'], row['BAD']
    return 1 / np.sqrt(calc_binom_variance(n, B).mean())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add annotations to the bed file.')
    parser.add_argument('bed_file', type=str, help='Path to the bed file')
    parser.add_argument('output', type=str, help='Path to the output file')
    parser.add_argument('--hotspots', type=str, help='Path to the in-hotspot indicator file')
    parser.add_argument('--footprints', type=str, help='Path to the in-footprint indicator file')
    parser.add_argument('--peaks', type=str, help='Path to the in-peaks indicator file')
    args = parser.parse_args()

    print('Reading files')

    df = pd.read_table(args.bed_file)
    initial_df_len = len(df.index)
    hotspots = read_indicator(args.hotspots)
    footprints = read_indicator(args.footprints)
    peaks = read_indicator(args.peaks)

    assert len(peaks) == len(hotspots) == len(footprints) == initial_df_len, f"Hotspots {len(hotspots)}, footprints {len(footprints)} and peaks {len(peaks)} should have the same length as the bed file {initial_df_len}."
    df['hotspots'] = hotspots
    df['footprints'] = footprints
    df['peaks'] = peaks
    
    data = df[['BAD', 'coverage']].drop_duplicates()
    data['inverse_mse'] = data.progress_apply(calc_inverse_mse, axis=1)
    df = df.drop(columns=['inverse_mse', 'inverse_mse_x', 'inverse_mse_y'], errors='ignore').merge(data, on=['BAD', 'coverage'])
    assert initial_df_len == len(df.index)
    df.to_csv(args.output, index=False, sep='\t')
    