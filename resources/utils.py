import subprocess
from os.path import *
import pybel


def inchi2smi(inchi, verbose=False):
    '''Converts InChI string to SMILES string.'''

    res = subprocess.check_output('echo "%s" | obabel -iinchi -ocan' % inchi,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')
    if verbose:
        print(res)

    res = [x.strip() for x in res.split('\n') if x is not '']

    if 'molecule converted' in res[-1]:
        return res[-2]

    return None


def smi2inchi(smi, verbose=False):
    '''Converts SMILES string to InChI string.'''

    res = subprocess.check_output('obabel -:"%s" -oinchi' % smi,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')
    if verbose:
        print(res)

    res = [x.strip() for x in res.split('\n') if x is not '']

    if 'molecule converted' in res[-1]:
        return res[-2]

    return None


def read_string(path):
    '''Reads first line from a file.'''

    with open(path, 'r') as f:
        return f.readlines()[0]


def write_string(string, path):
    '''Writes a string to file.'''

    with open(path, 'w') as f:
        f.write(string.strip() + '\n')


def desalt(inchi, verbose=False):
    '''Desalts an InChI string.'''

    smi = inchi2smi(inchi, verbose=verbose)

    if smi is None:
        return None

    smi = smi.replace('\"', '')

    return smi2inchi(smi, verbose=verbose)


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


def major_tautomer(inchi, verbose=False):
    '''Determines major tautomer of InChI string.'''

    res = subprocess.check_output('cxcalc majortautomer -f inchi "%s"' % inchi,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')
    if verbose:
        print(res)

    res = [x.strip() for x in res.split('\n') if x is not '']
    for line in res:
        if line.startswith('InChI='):
            return line
    return None


def inchi2formula(inchi, verbose=False):
    '''Determines formula from InChI string.'''

    res = subprocess.check_output('cxcalc formula "%s"' % inchi,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')
    if verbose:
        print(res)

    res = [x.strip() for x in res.split('\n') if x is not '']
    return res[-1].split()[-1].strip()


def inchi2geom(inchi, outfile, pngfile, ffield='gaff'):
    '''Converts InChI string to .mol geometry and saves a 2D visualization.'''

    mol = pybel.readstring("inchi", inchi)
    mol.addh()  # not necessary, because pybel make3D will add hydrogen

    mol.draw(show=False, filename=pngfile,
             usecoords=False, update=False)

    # Optimize 3D geometry of the molecule using pybel's make3D()
    mol.make3D(forcefield=ffield, steps=50)
    mol.localopt(forcefield=ffield, steps=500)

    mol.write('mol', outfile, True)


def calculatepKa(mol):
    '''Calculate pKa from .mol file.'''

    res = subprocess.check_output('cxcalc pka -i -40 -x 40 -d large %s' % mol,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')
    return res


def read_mass(path):
    '''Reads mass from molmass.py output file.'''

    with open(path, 'r') as f:
        lines = f.readlines()
        for x in lines:
            if 'Monoisotopic mass' in x:
                return float(x.split()[-1])


def read_pka(path):
    '''Reads pKa from cxcalc output'''

    with open(path, 'r') as f:
        vals = f.readlines()[1].split('\t')

    pk = [float(x) for x in vals[1:5]]
    atoms = [int(x) for x in vals[-1].split(',')]
    label = ['pka', 'pka', 'pkb', 'pkb']

    return [[a, p, x] for a, p, x in zip(atoms, pk, label)]
