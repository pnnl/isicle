import argparse
import pandas as pd
import numpy as np
from isicle import __version__


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform CCS calibration')
    parser.add_argument('infile', help='path to input .tsv file')
    parser.add_argument('outfile', help='path to output .tsv file')
    parser.add_argument('m', type=float, help='slope calibration')
    parser.add_argument('b', type=float, help='intercept calibration')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='Print version and exit.')

    args = parser.parse_args()
    args.adduct = args.adduct[1:-1]

    df = pd.read_csv(args.infile, sep='\t')
    df['ccs'] = args.m * df['ccs'] + args.b
    df['ccs_std'] = np.abs(args.m) * df['ccs_std']

    df.to_csv(args.outfile, sep='\t', index=False)
