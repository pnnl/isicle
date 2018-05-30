import subprocess
from os.path import *
import pybel


def inchi2smi(inchi):
    return subprocess.check_output('echo "%s" | obabel -iinchi -ocan' % inchi,
                                   stderr=subprocess.STDOUT, shell=True).decode('ascii').split('\n')[0].strip()


def smi2inchi(smi):
    return subprocess.check_output('obabel -:"%s" -oinchi' % smi,
                                   stderr=subprocess.STDOUT, shell=True).decode('ascii').split('\n')[0].strip()


def read_string(path):
    with open(path, 'r') as f:
        return f.readlines()[0]


def write_string(string, path):
    with open(path, 'w') as f:
        f.write(string + '\n')


def desalt(inchi):
    smi = inchi2smi(inchi)

    if smi is None:
        return None

    smi = smi.replace('\"', '')

    return smi2inchi(smi)


def neutralize(inchi):
    if 'q' in inchi:
        layers = inchi.split('/')
        new = layers[0]
        for i in range(1, len(layers)):
            if 'q' not in layers[i]:
                new += '/' + layers[i]
        return new
    return inchi


def major_tautomer(inchi):
    # timeout?
    return subprocess.check_output('cxcalc majortautomer -f inchi "%s"' % inchi,
                                   stderr=subprocess.STDOUT, shell=True).decode('ascii').split('\n')[0].strip()


def inchi2formula(inchi):
    return subprocess.check_output('cxcalc formula "%s"' % inchi,
                                   stderr=subprocess.STDOUT, shell=True).decode('ascii').split('\n')[1].split()[-1].strip()


def inchi2geom(inchi, mol2D, mol3D, xyz, ffield='gaff'):
    mol = pybel.readstring("inchi", inchi)
    mol.addh()  # not necessary, because pybel make3D will add hydrogen

    mol.write('mol', mol2D, True)
    mol.draw(show=False, filename=splitext(mol2D)[0] + '.png',
             usecoords=False, update=False)

    # Optimize 3D geometry of the molecule using pybel's make3D()
    mol.make3D(forcefield=ffield, steps=50)
    mol.localopt(forcefield=ffield, steps=500)

    mol.write('mol', mol3D, True)
    mol.write("xyz", xyz, True)


def read_mass(path):
    with open(path, 'r') as f:
        lines = f.readlines()
        for x in lines:
            if 'Monoisotopic mass' in x:
                return float(x.split()[-1])


if __name__ == '__main__':
    inchi = desalt('InChI=1S/C18H27N/c1-2-3-15-19(18-12-8-5-9-13-18)16-14-17-10-6-4-7-11-17/h2,4,6-7,10-11,18H,1,3,5,8-9,12-16H2')
    inchi = neutralize(inchi)
    inchi = major_tautomer(inchi)
    print(inchi)
    formula = inchi2formula(inchi)
    print(formula)
