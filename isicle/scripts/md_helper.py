import argparse
import numpy as np
from isicle.utils import pop_atom, push_atom
import shutil
from os.path import *
from isicle import __version__


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
    parser = argparse.ArgumentParser(description='Prepare molecule for MD')
    parser.add_argument('infile', help='path to .mol2 file')
    parser.add_argument('outfile', help='path to output .mol2 file')
    parser.add_argument('adduct', help='adduct type')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--prepare', action='store_true', help='prepare mode')
    mode.add_argument('--restore', action='store_true', help='restore mode')

    args = parser.parse_args()
    args.adduct = args.adduct[1:-1]

    if args.prepare is True:
        prepare(args.infile, args.outfile, adduct=args.adduct)
    elif args.restore is True:
        restore(args.infile, args.outfile, adduct=args.adduct)
