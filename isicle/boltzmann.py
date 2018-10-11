import argparse


__version__ = '0.1.0'


def boltzmann(infile, outfile):
    df = pd.read_csv(infile, sep='\t')

    g = df['DFT Energy'].values * 627.503
    mn = g.min()
    relG = g - mn
    b = np.exp(-relG / 0.5924847535)
    w = (b / b.sum()) * len(b)

    ws = DescrStatsW(df['Mean CCS'], weights=w, ddof=0)

    res = pd.Series([ws.mean, ws.std, ws.std_mean, ws.var, len(df.index)],
                    index=['mean', 'std', 'std_mean', 'var', 'N'])

    res.to_csv(outfile, sep='\t', index=False, header=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform Boltzmann weighting by energy across conformers.')
    parser.add_argument('infile', help='Path to input .tsv file.')
    parser.add_argument('outfile', help='Path to output .tsv file.')
    parser.add_argument('--version', action='version', version=__version__, help='Print version and exit.')

    args = parser.parse_args()

    import pandas as pd
    import numpy as np
    from statsmodels.stats.weightstats import DescrStatsW

    boltzmann(args.infile, args.outfile)
