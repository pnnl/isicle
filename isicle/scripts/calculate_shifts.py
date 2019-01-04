import argparse
import pandas as pd
from isicle import __version__


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate chemical shifts from shielding values')
    parser.add_argument('infile', help='path to shielding .tsv file')
    parser.add_argument('ref', help='path to reference .tsv file')
    parser.add_argument('outfile', help='path to output .tsv file')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    args = parser.parse_args()

    df = pd.read_csv(args.infile, sep='\t')
    ref = pd.read_csv(args.ref, sep='\t')

    df['shift'] = ref['shielding'].values - df['shielding'].values
    df['shift_std'] = (ref['shielding_std'].values ** 2 + df['shielding_std'].values ** 2) ** 0.5
    df = df['index', 'atom', 'shift', 'shift_std', 'shielding', 'shielding_std', 'n']

    df.to_csv(args.outfile, sep='\t', index=False)
