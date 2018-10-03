import subprocess
import pandas as pd
from os.path import *
import pybel
from . import geometry
import numpy as np
import platform
import openbabel


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

    try:
        res = subprocess.check_output(cmd,
                                      stderr=subprocess.STDOUT, shell=True).decode('ascii')
    except:
        return None

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

    return df['CCS_TJM'].values[0]


def read_mol(path, fmt='mol2'):
    return Mol(next(pybel.readfile(fmt, path)))


class Mol(pybel.Molecule):
    def __init__(self, mol):
        super().__init__(mol)

    def total_partial_charge(self):
        return np.array([a.partialcharge for a in self.atoms]).sum()

    def natoms(self):
        return len(self.atoms)


def pop_atom(path, output, atom='Na'):
    to_remove = []
    to_save = []
    with open(path, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if i == 2:
                info = [int(x) for x in line.split()]
            elif atom.upper() in line:
                to_remove.append(i)
                to_save.append(line)

    change = len(to_remove)
    with open(output, 'w') as f:
        for i, line in enumerate(lines):
            if i == 2:
                info[0] -= change
                f.write(' %s %s %s %s %s\n' % tuple(info))
            elif i in to_remove:
                pass
            else:
                f.write(line)

    return to_remove, to_save


def push_atom(path, output, idx, content):
    with open(path, 'r') as f:
        lines = f.readlines()

    info = [int(x) for x in lines[2].split()]

    for i, line in zip(idx, content):
        lines.insert(i, line)

    change = len(idx)
    with open(output, 'w') as f:
        for i, line in enumerate(lines):
            if i == 2:
                info[0] += change
                f.write(' %s %s %s %s %s\n' % tuple(info))
            else:
                f.write(line)


def select_frames(path, frames=10, low=1.25E6, high=1.45E6):
    tmp = []
    with open(path, 'r') as f:
        for line in f:
            if 'NSTEP' in line:
                props = [x for x in line.split() if x is not '=']
                d = {'step': int(props[1]),
                     'time': float(props[3]),
                     'temp': float(props[5])}
                tmp.append(d)
    df = pd.DataFrame(tmp)
    stepsize = df['step'][1] - df['step'][0]
    df['frame'] = df['step'] // stepsize

    ss = df[(df['step'] >= low) & (df['step'] <= high)]

    return ss['frame'].values[-frames:]


def standardizeMol2(path, reference, output):
    traj = next(pybel.readfile("mol2", path))
    ref = next(pybel.readfile("mol2", reference))
    for iatom in traj.atoms:
        ob = iatom.OBAtom
        idx = ob.GetIndex()

        jatom = ref.OBMol.GetAtomById(idx)
        ob.SetType(jatom.GetType())
        ob.SetAtomicNum(jatom.GetAtomicNum())

    traj.write("xyz", output, True)


def rmsd(mol1, mol2):
    a = next(pybel.readfile("xyz", mol1))
    b = next(pybel.readfile("xyz", mol2))

    align = openbabel.OBAlign(False, True)

    align.SetRefMol(a.OBMol)
    align.SetTargetMol(b.OBMol)
    align.Align()
    return align.GetRMSD()


def cycles(n):
    return ['%03d' % x for x in range(1, n + 1)]


def frames(n):
    return ['%03d' % x for x in range(n)]
