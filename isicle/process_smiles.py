import argparse
import subprocess
from core.utils import read_string, write_string
import sys


__version__ = '0.1.0'


def inchi2smi(inchi):
    '''Converts InChI string to canonical SMILES string.'''

    try:
        res = subprocess.check_output('echo "%s" | obabel -iinchi -ocan' % inchi,
                                      stderr=subprocess.STDOUT, shell=True).decode('ascii')
    except:
        return None

    res = [x.strip() for x in res.split('\n') if x is not '']

    if 'molecule converted' in res[-1]:
        return res[-2]

    return None


def canonicalize(smiles):
    '''Converts SMILES string to canonical SMILES string.'''

    try:
        res = subprocess.check_output('echo "%s" | obabel -ismi -ocan' % smiles,
                                      stderr=subprocess.STDOUT, shell=True).decode('ascii')
    except:
        return None

    res = [x.strip() for x in res.split('\n') if x is not '']

    if 'molecule converted' in res[-1]:
        return res[-2]

    return None


def desalt(smiles):
    '''Desalts a canonical SMILES string.'''

    try:
        res = subprocess.check_output('echo "%s" | obabel -ican -r -ocan' % smiles,
                                      stderr=subprocess.STDOUT, shell=True).decode('ascii')
    except:
        return None

    res = [x.strip() for x in res.split('\n') if x is not '']

    if 'molecule converted' in res[-1]:
        return res[-2]

    return None


def neutralize(smiles):
    '''Neutralizes an canonical SMILES string.'''

    # currently not operational for SMILES
    return smiles


def tautomerize(smiles):
    '''Determines major tautomer of canonical SMILES string.'''

    res = subprocess.check_output('cxcalc majortautomer -f smiles "%s"' % smiles,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')

    res = [x.strip() for x in res.split('\n') if x is not '']

    if len(res) > 1:
        return None
    else:
        return canonicalize(res[0])


def smiles2formula(smiles):
    '''Determines formula from canonical SMILES string.'''

    res = subprocess.check_output('cxcalc formula "%s"' % smiles,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')

    res = [x.strip() for x in res.split('\n') if x is not '']
    return res[-1].split()[-1].strip()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process InChI/SMILES string.')
    parser.add_argument('infile', help='Path to input InChI (.inchi) or canonical SMILES (.smi) file.')
    parser.add_argument('outfile', help='Path to output canonical SMILES (.smi) file.')

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--inchi', action='store_true', help='InChI to SMILES mode.')
    mode.add_argument('--desalt', action='store_true', help='Desalt mode.')
    mode.add_argument('--neutralize', action='store_true', help='Neutralize mode.')
    mode.add_argument('--tautomerize', action='store_true', help='Tautomerize mode.')
    mode.add_argument('--formula', action='store_true', help='Formula mode.')

    args = parser.parse_args()

    s = read_string(args.infile)

    if args.inchi is True:
        smiles = inchi2smi(s)
    elif args.desalt is True:
        smiles = desalt(s)
    elif args.neutralize is True:
        smiles = neutralize(s)
    elif args.tautomerize is True:
        smiles = tautomerize(s)
    elif args.formula is True:
        smiles = smiles2formula(s)

    if smiles is not None:
        write_string(smiles, args.outfile)
    else:
        sys.exit(1)
