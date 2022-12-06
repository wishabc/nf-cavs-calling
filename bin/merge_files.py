import shutil
import sys

paths = sys.argv[2:]

with open(sys.argv[1], 'wb') as outfile:
    for i, fname in enumerate(paths):
        with open(fname, 'rb') as infile:
            if i != 0:
                infile.readline()  # Throw away header on all but first file
            # Block copy rest of file from input to output without parsing
            shutil.copyfileobj(infile, outfile)
