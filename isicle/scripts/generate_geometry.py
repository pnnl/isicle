import argparse
from openbabel import pybel
from isicle.utils import read_string
from os.path import *
from isicle import __version__


def inchi2geom(inchi, forcefield='mmff94', steps=500):
    '''Converts InChI string to .mol geometry and saves a 2D visualization.'''

    mol = pybel.readstring("inchi", inchi)

    # Optimize 3D geometry of the molecule using pybel's make3D()
    mol.make3D(forcefield=forcefield, steps=50)
    mol.localopt(forcefield=forcefield, steps=steps)

    return mol


def smiles2geom(smiles, forcefield='mmff94', steps=500):
    '''Converts canonical SMILES string to .mol geometry and saves a 2D visualization.'''

    mol = pybel.readstring("can", smiles)

    # Optimize 3D geometry of the molecule using pybel's make3D()
    mol.make3D(forcefield=forcefield, steps=50)
    mol.localopt(forcefield=forcefield, steps=steps)

    return mol


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate geometry from InChI string')
    parser.add_argument('infile', help='path to input inchi (.inchi) or smiles (.smi) file')
    parser.add_argument('mol', help='path to output .mol file')
    parser.add_argument('mol2', help='path to output .mol2 file')
    parser.add_argument('xyz', help='path to output .xyz file')
    parser.add_argument('png', help='path to output .png image')
    parser.add_argument('-ff', '--forcefield', type=str, default='gaff', help='forcefield type')
    parser.add_argument('-s', '--steps', type=int, default=500, help='number of forcefield optimization steps')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    args = parser.parse_args()

    s = read_string(args.infile)

    if splitext(args.infile)[-1].lower() == '.inchi':
        mol = inchi2geom(s, forcefield=args.forcefield,
                         steps=args.steps)
    elif splitext(args.infile)[-1].lower() == '.smi':
        mol = smiles2geom(s, forcefield=args.forcefield,
                          steps=args.steps)
    else:
        raise IOError('File type "%s" not supported.' % splitext(args.infile)[-1].lower())

    mol.draw(show=False, filename=args.png, usecoords=False, update=False)
    mol.write('mol', args.mol, overwrite=True)
    mol.write('mol2', args.mol2, overwrite=True)
    mol.write('xyz', args.xyz, overwrite=True)
