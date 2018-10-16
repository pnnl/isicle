import argparse


__version__ = '0.1.0'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate chemical shifts from shielding values.')
    parser.add_argument('infile', help='Paths to shielding .tsv file.')
    parser.add_argument('ref', help='Path to reference .tsv file.')
    parser.add_argument('outfile', help='Path to output .tsv file.')
    parser.add_argument('--version', '-v', action='version', version=__version__, help='Print version and exit.')

    args = parser.parse_args()

    import pandas as pd

    df = pd.read_tsv(args.infile, sep='\t')
    ref = pd.read_tsv(args.ref, sep='\t')

    df['mean'] = ref['mean'].values - df['mean'].values
    df['std'] = (ref['std'].values ** 2 + df['std'].values ** 2) ** 0.5
    df = df['index', 'atom', 'mean', 'std', 'N']

    df.to_csv(args.outfile, sep='\t', index=False)
