import subprocess
import pandas as pd
import pybel
import numpy as np
import platform


def tail(f, lines=1):
    return [x.strip() for x in subprocess.check_output(['tail', '-n%s' % lines, f]).decode('ascii').split('\n')]


def read_string(path):
    '''Reads first line from a file.'''

    with open(path, 'r') as f:
        return f.readlines()[0].strip()


def write_string(string, path):
    '''Writes a string to file.'''

    with open(path, 'w') as f:
        f.write(string.strip() + '\n')


def read_mass(path):
    '''Reads mass from molmass.py output file.'''

    with open(path, 'r') as f:
        lines = f.readlines()
        for x in lines:
            if 'Monoisotopic mass' in x:
                return float(x.split()[-1])


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
        if 'NA' in line:
            parts = line.split()
            parts[1] = 'NA'
            parts[5] = 'Na+'
            parts[6] = '2'
            parts[7] = 'Na+'
            parts[8] = '1.0000'
            line = ' '.join(parts) + '\n'
        lines.insert(i + 1, line)

    change = len(idx)
    with open(output, 'w') as f:
        for i, line in enumerate(lines):
            if i == 2:
                info[0] += change
                f.write(' %s %s %s %s %s\n' % tuple(info))
            else:
                f.write(line)


def getOS():
    system = platform.system().lower()
    if system == 'darwin':
        return 'osx'
    return system


def cycles(n):
    return ['%03d' % x for x in range(1, n + 1)]


def frames(n):
    return ['%03d' % x for x in range(n)]
