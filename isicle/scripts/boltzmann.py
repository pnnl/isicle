import argparse
import pandas as pd
import numpy as np
from statsmodels.stats.weightstats import DescrStatsW
from isicle import __version__


def ccs(infile, outfile):
    df = pd.read_csv(infile, sep='\t')

    g = df['dft_energy'].values * 627.503
    mn = g.min()
    relG = g - mn
    b = np.exp(-relG / 0.5924847535)
    w = (b / b.sum()) * len(b)

    ws = DescrStatsW(df['ccs'], weights=w, ddof=0)

    res = pd.DataFrame([[ws.mean, ws.std, len(df.index)]],
                       columns=['ccs', 'ccs_std', 'n'])

    res.to_csv(outfile, sep='\t', index=False, header=True)


def shielding(infile, outfile):
    df = pd.read_csv(infile, sep='\t')

    data = []
    for name, group in df.groupby(['index', 'atom']):
        g = group['dft_energy'].values * 627.503
        mn = g.min()
        relG = g - mn
        b = np.exp(-relG / 0.5924847535)
        w = (b / b.sum()) * len(b)

        ws = DescrStatsW(group['shielding'], weights=w, ddof=0)
        data.append([name[0], name[1], ws.mean, ws.std, len(group.index)])

    df2 = pd.DataFrame(data, columns=['index', 'atom', 'shielding', 'shielding_std', 'n'])
    df2.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform Boltzmann weighting by energy across conformers')
    parser.add_argument('infile', help='path to input .tsv file')
    parser.add_argument('outfile', help='path to output .tsv file')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--ccs', action='store_true', help='ccs mode')
    mode.add_argument('--shielding', action='store_true', help='shielding mode')

    args = parser.parse_args()

    if args.ccs is True:
        ccs(args.infile, args.outfile)
    elif args.shielding is True:
        shielding(args.infile, args.outfile)
