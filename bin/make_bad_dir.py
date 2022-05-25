from helpers import make_bad_dir
import shutil
import argparse
import os

def main(bad, stats_file_path):
    dir = os.path.dirname(stats_file_path)
    _, bad_dir = make_bad_dir(dir, bad)
    shutil.move(stats_file_path, bad_dir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create dir for negbinfit')
    parser.add_argument('-s', help='stats.tsv file path')
    parser.add_argument('-b', help='BAD score')
    args = parser.parse_args()
    main(args.b, args.s)