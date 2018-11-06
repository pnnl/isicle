import argparse
import pandas as pd
from isicle.utils import tail
from isicle import __version__


def parse_mobcal(f):
    done = False

    lines = tail(f, lines=10)

    for line in lines:
        if "average (second order) TM mobility" in line:
            m_mn = float(line.split('=')[-1])
        elif "average TM cross section" in line:
            ccs_mn = float(line.split('=')[-1])
        elif "standard deviation" in line:
            ccs_std = float(line.split('=')[-1])
            done = True
    if done is True:
        return [m_mn, ccs_mn, ccs_std]
    else:
        return None


def batch(ccsfiles, efiles, output):
    res = []
    for ccsfile, efile in zip(ccsfiles, efiles):
        tmp = parse_mobcal(ccsfile)
        if tmp is not None:
            # get energy
            with open(efile, 'r') as f:
                e = float(f.readlines()[0])
            tmp.append(e)
            res.append(tmp)

    df = pd.DataFrame(res, columns=['mobility', 'ccs', 'ccs_std', 'dft_energy'])
    df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse mobcal output')
    parser.add_argument('outfile', help='path to output file')
    parser.add_argument('--infiles', nargs='+', required=True, help='paths to mobcal .out files')
    parser.add_argument('--efiles', nargs='+', required=True, help='paths to .energy files')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    args = parser.parse_args()

    assert len(args.infiles) == len(args.efiles), 'Number of output and energy files must be equal.'

    args.infiles.sort()
    args.efiles.sort()

    batch(args.infiles, args.efiles, args.outfile)
