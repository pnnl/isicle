import argparse
import subprocess
from os.path import *


__version__ = '0.1.0'


def smi2key(smi):
    try:
        res = subprocess.check_output('obabel -:"%s" -oinchikey' % smi,
                                      stderr=subprocess.STDOUT, shell=True).decode('ascii')
    except:
        return None

    res = [x.strip() for x in res.split('\n') if x is not '']

    if 'molecule converted' in res[-1]:
        return res[-2]

    return None


def inchi2key(inchi):
    try:
        res = subprocess.check_output('echo "%s" | obabel -iinchi -oinchikey' % inchi,
                                      stderr=subprocess.STDOUT, shell=True).decode('ascii')
    except:
        return None

    res = [x.strip() for x in res.split('\n') if x is not '']

    if 'molecule converted' in res[-1]:
        return res[-2]

    return None


def cli():
    parser = argparse.ArgumentParser(description="A program to convert a .txt file of InChI or SMILES to individual files for ISiCLE.")
    parser.add_argument('filepath', type=argparse.FileType('r'), help="Include valid filepath to .txt file.")
    parser.add_argument('--config', required=True, help="Path to ISiCLE configuration file.")
    parser.add_argument('--version', '-v', action='version', version=__version__, help='Print version and exit.')

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--smi', action='store_true', help='SMILES mode.')
    mode.add_argument('--inchi', action='store_true', help='InChI mode.')

    args = parser.parse_args()

    import pandas as pd
    import yaml
    from isicle.utils import write_string
    import os

    df = pd.read_csv(args.filepath, sep='\n', header=None)

    with open(args.config, 'r') as f:
        outdir = join(yaml.load(f)['path'], 'input')

    if not exists(outdir):
        os.makedirs(outdir)

    if args.smi is True:
        for row in df.values:
            smi = row[0]
            key = smi2key(smi)
            if key is not None:
                write_string(smi, join(outdir, '%s.smi') % key)
            else:
                print('Failed to hash %s.' % smi)

    elif args.inchi is True:
        for row in df.values:
            inchi = row[0]
            key = inchi2key(inchi)
            if key is not None:
                write_string(inchi, join(outdir, '%s.inchi') % key)
            else:
                print('Failed to hash %s.' % inchi)


if __name__ == '__main__':
    cli()
