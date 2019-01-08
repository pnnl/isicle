import argparse
import pandas as pd
import numpy as np
from isicle import __version__


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate chemical shifts from shielding values')
    parser.add_argument('infile', help='path to shielding .tsv file')
    parser.add_argument('ref', help='reference shielding dictionary')
    parser.add_argument('outfile', help='path to output .tsv file')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    args = parser.parse_args()

    df = pd.read_csv(args.infile, sep='\t')

    df['shift'] = np.nan
    df['shift_std'] = np.nan

    for k, v in eval(args.ref).items():
        df.loc[df['atom'] == k, 'shift'] = v['mean'] - df.loc[df['atom'] == k, 'shielding']
        df.loc[df['atom'] == k, 'shift_std'] = np.sqrt(np.square(df.loc[df['atom'] == k, 'shielding_std']) + np.square(v['std']))

    df = df[['index', 'atom', 'shift', 'shift_std', 'n']]

    df.to_csv(args.outfile, sep='\t', index=False)
