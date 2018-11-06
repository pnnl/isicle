import argparse
import pandas as pd
from isicle import __version__


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine shielding data')
    parser.add_argument('infiles', nargs='+', help='paths to .shielding files')
    parser.add_argument('outfile', help='path to output .tsv file')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    args = parser.parse_args()

    dfs = [pd.read_csv(x, sep='\t') for x in args.infiles]

    df = pd.concat(dfs, axis=0, ignore_index=True)

    df.to_csv(args.outfile, sep='\t', index=False)
