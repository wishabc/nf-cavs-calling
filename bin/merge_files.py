import shutil
import sys

paths = sys.argv[3:]
group_key = sys.argv[2]

with open(sys.argv[1], 'w') as outfile:
    for i, fname in enumerate(paths):
        with open(fname, 'r') as infile:
            for j, line in enumerate(infile):
                if i == 0 and j == 0:
                    new_line = line.replace('\n', '\tgroup_id\n')
                elif i != 0 and j == 0:
                    continue
                else:
                    new_line = line.replace('\n', f'\t{group_key}\n')
                outfile.write(new_line)
