import argparse
from core.utils import read_string
import numpy as np
from os.path import *
import shutil

__version__ = '0.1.0'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Downselect conformers by RMSD.')
    parser.add_argument('infiles', nargs='+', help='Paths to .xyz files.')
    parser.add_argument('rfiles', nargs='+', help='Paths to .rmsd files.')
    parser.add_argument('outdir', help='Path to output directory.')

    args = parser.parse_args()

    assert len(args.infiles) == len(args.rfiles), 'Number of .xyz and .rmsd files must be equal.'

    args.infiles.sort()
    args.rfiles.sort()

    vals = []
    for f in args.rfiles:
        vals.append(float(read_string(f)))

    vals = np.array(vals)
    idx = np.argsort(vals)

    s = args.infiles[idx[0]]
    d1 = args.infiles[idx[-1]]
    d2 = args.infiles[idx[-2]]

    sout = join(outdir, basename(s).rsplit('_', 1)[0] + '_s.xyz')
    d1out = join(outdir, basename(d1).rsplit('_', 1)[0] + '_d1.xyz')
    d2out = join(outdir, basename(d2).rsplit('_', 1)[0] + '_d2.xyz')

    shutil.copy2(s, sout)
    shutil.copy2(d1, d1out)
    shutil.copy2(d2, d2out)
