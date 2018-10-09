import argparse
import pybel
from core.utils import read_string

__version__ = '0.1.0'


def inchi2geom(inchi, forcefield='mmff94', steps=500):
    '''Converts InChI string to .mol geometry and saves a 2D visualization.'''

    mol = pybel.readstring("inchi", inchi)
    mol.addh()  # not necessary, because pybel make3D will add hydrogen

    # Optimize 3D geometry of the molecule using pybel's make3D()
    mol.make3D(forcefield=forcefield, steps=50)
    mol.localopt(forcefield=forcefield, steps=steps)

    return mol


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate geometry from InChI string.')
    parser.add_argument('infile', help='Path to input .inchi file.')
    parser.add_argument('mol', help='Path to output .mol file.')
    parser.add_argument('mol2', help='Path to output .mol2 file.')
    parser.add_argument('xyz', help='Path to output .xyz file.')
    parser.add_argument('png', help='Path to output .png image.')
    parser.add_argument('--forcefield', '-ff', type=str, default='gaff', help='Forcefield type.')
    parser.add_argument('--steps', '-s', type=int, default=500, help='Number of forcefield optimization steps.')

    args = parser.parse_args()

    inchi = read_string(args.infile)
    mol = inchi2geom(inchi, forcefield=args.forcefield,
                     steps=args.iterations)

    mol.draw(show=False, filename=args.png, usecoords=False, update=False)
    mol.write('mol', args.mol, overwrite=True)
    mol.write('xyz', args.xyz, overwrite=True)
    mol.write('mol2', args.mol2, overwrite=True)
