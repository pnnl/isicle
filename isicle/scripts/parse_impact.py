import argparse
import pandas as pd
from isicle.utils import read_mass, write_string
from isicle import __version__


def read_impact(path):
    # read ccs file
    df = pd.read_csv(path, delim_whitespace=True, index_col=False)

    # clean up
    df.drop(['#Str', 'nr', 'filename', '(SEM_rel)'], axis=1, inplace=True)
    df.columns = ['CCS_PA', 'SEM_rel', 'CCS_TJM']
    df['SEM_rel'] = float(df['SEM_rel'].str[:-1])

    return df['CCS_TJM'].values[0]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse IMPACT output')
    parser.add_argument('infile', help='path to impact .out files')
    parser.add_argument('mass', help='path to mass files.')
    parser.add_argument('out_He', help='path to helium ccs output file')
    parser.add_argument('out_N2', help='path to nitrogen ccs output file')
    parser.add_argument('-a', '--alpha', type=float, default=27.9, help='alpha calibration')
    parser.add_argument('-b', '--beta', type=float, default=0.14, help='beta calibration')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    args = parser.parse_args()

    # read inputs
    ccs_He = read_impact(args.infile)
    m = read_mass(args.mass)

    ccs_N2 = ccs_He + args.alpha * m ** args.beta

    write_string(str(ccs_He), args.out_He)
    write_string(str(ccs_N2), args.out_N2)
