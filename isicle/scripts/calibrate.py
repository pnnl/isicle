import argparse


__version__ = '0.1.0'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform CCS calibration.')
    parser.add_argument('infile', help='Path to input .tsv file.')
    parser.add_argument('outfile', help='Path to output .tsv file.')
    parser.add_argument('m', type=float, help='Slope calibration.')
    parser.add_argument('b', type=float, help='Intercept calibration.')
    parser.add_argument('--version', '-v', action='version', version=__version__, help='Print version and exit.')

    args = parser.parse_args()

    import pandas as pd
    import numpy as np

    args.adduct = args.adduct[1:-1]

    df = pd.read_csv(args.infile, sep='\t')
    df['ccs'] = args.m * df['ccs'] + args.b
    df['ccs_std'] = np.abs(args.m) * df['ccs_std']

    df.to_csv(args.outfile, sep='\t', index=False)
