import subprocess
from os.path import *
import pybel
from . import geometry
import pandas as pd
import numpy as np
import platform
import sys


def getOS():
    system = platform.system().lower()
    if system == 'darwin':
        return 'osx'
    return system


def inchi2smi(inchi, desalt=False, log=None):
    '''Converts InChI string to SMILES string.'''

    if desalt is True:
        ds = '-r'
        sep = '\n'
    else:
        ds = ''
        sep = ''

    cmd = 'echo "%s" | obabel -iinchi %s -ocan' % (inchi, ds)
    res = subprocess.check_output(cmd,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')
    if log is not None:
        with open(log, 'a') as f:
            f.write('inchi2smi:\n')
            f.write(cmd + '\n\n')
            f.write(res + sep)

    res = [x.strip() for x in res.split('\n') if x is not '']

    if 'molecule converted' in res[-1]:
        return res[-2]

    return None


def smi2inchi(smi, log=None):
    '''Converts SMILES string to InChI string.'''

    cmd = 'obabel -:"%s" -oinchi' % smi
    res = subprocess.check_output(cmd,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')
    if log is not None:
        with open(log, 'a') as f:
            f.write('smi2inchi:\n')
            f.write(cmd + '\n\n')
            f.write(res)

    res = [x.strip() for x in res.split('\n') if x is not '']

    if 'molecule converted' in res[-1]:
        return res[-2]

    return None


def read_string(path):
    '''Reads first line from a file.'''

    with open(path, 'r') as f:
        return f.readlines()[0].strip()


def write_string(string, path):
    '''Writes a string to file.'''

    with open(path, 'w') as f:
        f.write(string.strip() + '\n')


def desalt(inchi, log=None):
    '''Desalts an InChI string.'''
    if log is not None:
        with open(log, 'w') as f:
            f.write('desalt:\n\n')

    smi = inchi2smi(inchi, desalt=True, log=log)

    if smi is None:
        return None

    return smi2inchi(smi, log=log)


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


def tautomerize(inchi, log=None):
    '''Determines major tautomer of InChI string.'''

    cmd = 'cxcalc majortautomer -f inchi "%s"' % inchi
    res = subprocess.check_output(cmd,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')
    if log is not None:
        with open(log, 'w') as f:
            f.write('tautomerize:\n')
            f.write(cmd + '\n\n')
            f.write(res)

    res = [x.strip() for x in res.split('\n') if x is not '']
    for line in res:
        if line.startswith('InChI='):
            return line
    return None


def inchi2formula(inchi, log=None):
    '''Determines formula from InChI string.'''

    cmd = 'cxcalc formula "%s"' % inchi
    res = subprocess.check_output(cmd,
                                  stderr=subprocess.STDOUT, shell=True).decode('ascii')
    if log is not None:
        with open(log, 'w') as f:
            f.write('inchi2formula:\n')
            f.write(cmd + '\n\n')
            f.write(res)

    res = [x.strip() for x in res.split('\n') if x is not '']
    return res[-1].split()[-1].strip()


def inchi2geom(inchi, forcefield='mmff94', steps=500):
    '''Converts InChI string to .mol geometry and saves a 2D visualization.'''

    mol = pybel.readstring("inchi", inchi)
    mol.addh()  # not necessary, because pybel make3D will add hydrogen

    # Optimize 3D geometry of the molecule using pybel's make3D()
    mol.make3D(forcefield=forcefield, steps=50)
    mol.localopt(forcefield=forcefield, steps=steps)

    return mol


def read_mass(path):
    '''Reads mass from molmass.py output file.'''

    with open(path, 'r') as f:
        lines = f.readlines()
        for x in lines:
            if 'Monoisotopic mass' in x:
                return float(x.split()[-1])


def create_adduct(mol, adduct, idx, forcefield='mmff94', steps=500):
    if '-' in adduct:
        hidx = geometry.nearestHydrogen(mol, idx)
        adduct = geometry.removeAtomFromMol(mol, hidx)
    elif '+' in adduct:
        atom = adduct.split('+')[-1]
        if atom.lower() == 'na':
            adduct = geometry.addAtomToMol(mol, atom, idx, covalent=False)
        else:
            adduct = geometry.addAtomToMol(mol, atom, idx, covalent=True)

    # talk to Jamie about this:
    adduct.localopt(forcefield=forcefield, steps=steps)
    return adduct


def read_pka(path):
    '''Reads pKa from cxcalc output'''

    df = pd.read_csv(path, sep='\t')

    # offset indices because cxcalc is 1-based
    idx = [int(x) - 1 for x in df['atoms'].values[0].split(',')]
    pk = df.values[0][1:5]
    label = ['a1', 'a2', 'b1', 'b2']

    res = {}
    i = 0
    for x, p in zip(label, pk):
        if not np.isnan(p):
            res[x] = idx[i]
            i += 1

    return res


def read_impact(path):
    # read ccs file
    df = pd.read_csv(path, delim_whitespace=True, index_col=False)

    # clean up
    df.drop(['#Str', 'nr', 'filename', '(SEM_rel)'], axis=1, inplace=True)
    df.columns = ['CCS_PA', 'SEM_rel', 'CCS_TJM']
    df['SEM_rel'] = float(df['SEM_rel'].str[:-1])

    return df
