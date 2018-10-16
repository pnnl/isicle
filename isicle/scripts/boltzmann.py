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


def shielding(infile, outfile):
    df = pd.read_csv(infile, sep='\t')

    data = []
    for name, group in df.groupby(['Index', 'Atom']):
        g = group['DFT Energy'].values * 627.503
        mn = g.min()
        relG = g - mn
        b = np.exp(-relG / 0.5924847535)
        w = (b / b.sum()) * len(b)

        ws = DescrStatsW(group['Shielding'], weights=w, ddof=0)
        data.append([name[0], name[1], ws.mean, ws.std, ws.std_mean, ws.var, len(group.index)])

    df2 = pd.DataFrame(data, columns=['index', 'atom', 'mean', 'std', 'std_mean', 'var', 'N'])
    df2.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform Boltzmann weighting by energy across conformers.')
    parser.add_argument('infile', help='Path to input .tsv file.')
    parser.add_argument('outfile', help='Path to output .tsv file.')
    parser.add_argument('--version', '-v', action='version', version=__version__, help='Print version and exit.')

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--ccs', action='store_true', help='CCS mode.')
    mode.add_argument('--shielding', action='store_true', help='Shielding mode.')

    args = parser.parse_args()

    import pandas as pd
    import numpy as np
    from statsmodels.stats.weightstats import DescrStatsW

    if args.ccs is True:
        boltzmann(args.infile, args.outfile)
    elif args.shielding is True:
        shielding(args.infile, args.outfile)
