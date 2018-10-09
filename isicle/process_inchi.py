import argparse
import subprocess
from core.utils import read_string, write_string
import sys

__version__ = '0.1.0'


def inchi2smi(inchi, desalt=False):
    '''Converts InChI string to SMILES string.'''

    if desalt is True:
        ds = '-r'
    else:
        ds = ''

    try:
        res = subprocess.check_output('echo "%s" | obabel -iinchi %s -ocan' % (inchi, ds),
                                      stderr=subprocess.STDOUT, shell=True).decode('ascii')
    except:
        return None

    res = [x.strip() for x in res.split('\n') if x is not '']

    if 'molecule converted' in res[-1]:
        return res[-2]

    return None


def smi2inchi(smi):
    '''Converts SMILES string to InChI string.'''

    res = subprocess.check_output('obabel -:"%s" -oinchi' % smi,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')

    res = [x.strip() for x in res.split('\n') if x is not '']

    if 'molecule converted' in res[-1]:
        return res[-2]

    return None


def desalt(inchi):
    '''Desalts an InChI string.'''

    smi = inchi2smi(inchi, desalt=True)

    if smi is None:
        return None

    return smi2inchi(smi)


def neutralize(inchi):
    '''Neutralizes an InChI string.'''

    if 'q' in inchi:
        layers = inchi.split('/')
        new = layers[0]
        for i in range(1, len(layers)):
            if 'q' not in layers[i]:
                new += '/' + layers[i]
        return new
    return inchi


def tautomerize(inchi):
    '''Determines major tautomer of InChI string.'''

    res = subprocess.check_output('cxcalc majortautomer -f inchi "%s"' % inchi,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')

    res = [x.strip() for x in res.split('\n') if x is not '']
    for line in res:
        if line.startswith('InChI='):
            return line
    return None


def inchi2formula(inchi):
    '''Determines formula from InChI string.'''

    res = subprocess.check_output('cxcalc formula "%s"' % inchi,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')

    res = [x.strip() for x in res.split('\n') if x is not '']
    return res[-1].split()[-1].strip()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process InChI string.')
    parser.add_argument('infile', help='Path to input .inchi file.')
    parser.add_argument('outfile', help='Path to output .inchi file.')

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--desalt', action='store_true', help='Desalt mode.')
    mode.add_argument('--neutralize', action='store_true', help='Neutralize mode.')
    mode.add_argument('--tautomerize', action='store_true', help='Tautomerize mode.')
    mode.add_argument('--formula', action='store_true', help='Formula mode.')

    args = parser.parse_args()

    inchi = read_string(args.infile)

    if args.desalt is True:
        inchi = desalt(inchi)
    elif args.neutralize is True:
        inchi = neutralize(inchi)
    elif args.tautomerize is True:
        inchi = tautomerize(inchi)
    elif args.formula is True:
        inchi = inchi2formula(inchi)

    if inchi is not None:
        write_string(inchi, args.outfile)
    else:
        sys.exit(1)
