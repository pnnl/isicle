import argparse


__version__ = '0.1.0'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Combine chemical shift data.')
    parser.add_argument('infiles', nargs='+', help='Paths to .shielding files.')
    parser.add_argument('outfile', help='Path to output .tsv file.')
    parser.add_argument('--version', '-v', action='version', version=__version__, help='Print version and exit.')

    args = parser.parse_args()

    import pandas as pd

    dfs = [pd.read_csv(x, sep='\t') for x in args.infiles]

    df = pd.concat(dfs, axis=0, ignore_index=True)

    df.to_csv(args.outfile, sep='\t', index=False)
