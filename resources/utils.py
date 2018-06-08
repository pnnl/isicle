import subprocess
from os.path import *
import pybel


def inchi2smi(inchi):
    '''Converts InChI string to SMILES string.'''

    return subprocess.check_output('echo "%s" | obabel -iinchi -ocan' % inchi,
                                   stderr=subprocess.STDOUT, shell=True).decode('ascii').split('\n')[0].strip()


def smi2inchi(smi):
    '''Converts SMILES string to InChI string.'''

    return subprocess.check_output('obabel -:"%s" -oinchi' % smi,
                                   stderr=subprocess.STDOUT, shell=True).decode('ascii').split('\n')[0].strip()


def read_string(path):
    '''Reads first line from a file.'''

    with open(path, 'r') as f:
        return f.readlines()[0]


def write_string(string, path):
    '''Writes a string to file.'''

    with open(path, 'w') as f:
        f.write(string + '\n')


def desalt(inchi):
    '''Desalts an InChI string.'''

    smi = inchi2smi(inchi)

    if smi is None:
        return None

    smi = smi.replace('\"', '')

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


def major_tautomer(inchi):
    '''Determines major tautomer of InChI string.'''

    return subprocess.check_output('cxcalc majortautomer -f inchi "%s"' % inchi,
                                   stderr=subprocess.STDOUT, shell=True).decode('ascii').split('\n')[0].strip()


def inchi2formula(inchi):
    '''Determines formula from InChI string.'''

    return subprocess.check_output('cxcalc formula "%s"' % inchi,
                                   stderr=subprocess.STDOUT, shell=True).decode('ascii').split('\n')[1].split()[-1].strip()


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


def read_mass(path):
    '''Reads mass from molmass.py output file.'''

    with open(path, 'r') as f:
        lines = f.readlines()
        for x in lines:
            if 'Monoisotopic mass' in x:
                return float(x.split()[-1])
