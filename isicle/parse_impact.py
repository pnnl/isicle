import argparse


__version__ = '0.1.0'


def read_impact(path):
    # read ccs file
    df = pd.read_csv(path, delim_whitespace=True, index_col=False)

    # clean up
    df.drop(['#Str', 'nr', 'filename', '(SEM_rel)'], axis=1, inplace=True)
    df.columns = ['CCS_PA', 'SEM_rel', 'CCS_TJM']
    df['SEM_rel'] = float(df['SEM_rel'].str[:-1])

    return df['CCS_TJM'].values[0]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse IMPACT output.')
    parser.add_argument('infile', help='Path to IMPACT .out files.')
    parser.add_argument('mass', help='Path to mass files.')
    parser.add_argument('out_He', help='Path to He CCS output file.')
    parser.add_argument('out_N2', help='Path to N2 CCS output file.')
    parser.add_argument('--alpha', '-a', type=float, default=27.9, help='Alpha calibration.')
    parser.add_argument('--beta', '-b', type=float, default=0.14, help='Beta calibration.')
    parser.add_argument('--version', '-v', action='version', version=__version__, help='Print version and exit.')

    args = parser.parse_args()

    import pandas as pd
    from isicle.core.utils import read_mass, write_string

    # read inputs
    ccs_He = read_impact(args.infile)
    m = read_mass(args.mass)

    ccs_N2 = ccs_He + args.alpha * m ** args.beta

    write_string(str(ccs_He), args.out_He)
    write_string(str(ccs_N2), args.out_N2)
