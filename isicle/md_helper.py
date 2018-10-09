import numpy as np
from core.utils import pop_atom
import shutil
import argparse
from os.path import *

__version__ = '0.1.0'


def prepare(infile, outfile, adduct='+H'):
    if adduct == '+Na':
        idx, content = pop_atom(infile, outfile, atom='Na')
    else:
        shutil.copy2(infile, outfile)
        idx = None
        content = None

    np.save(splitext(splitext(outfile)[0])[0] + '.idx.npy', idx)
    np.save(splitext(splitext(outfile)[0])[0] + '.content.npy', content)


def restore(infile, outfile, adduct='+H'):
    if adduct == '+Na':
        idx = np.load(splitext(splitext(infile)[0])[0] + '.idx.npy')
        content = np.load(splitext(splitext(infile)[0])[0] + '.content.npy')
        push_atom(infile, outfile, idx, content)
    else:
        shutil.copy2(infile, outfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare molecule for MD.')
    parser.add_argument('infile', help='Path to .mol2 file.')
    parser.add_argument('outfile', help='Path to output .mol2 file.')
    parser.add_argument('adduct', help='Adduct type.')

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true', help='Prepare mode.')
    mode.add_argument('--restore', action='store_true', help='Restore mode.')

    args = parser.parse_args()
    args.adduct = args.adduct[1:-1]

    if args.prepare is True:
        prepare(args.infile, args.outfile, adduct=args.adduct)
    elif args.restore is True:
        restore(args.infile, args.outfile, adduct=args.adduct)
