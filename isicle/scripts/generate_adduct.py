import argparse
import openbabel as ob
from isicle.resources import geometry
from isicle.utils import read_pka, read_mol, write_string
from isicle import __version__
# import subprocess
# from os.path import *


def create_adduct(mol, adduct, idx, forcefield='mmff94', steps=500):
    if adduct == '-H':
        hidx = geometry.nearestHydrogen(mol, idx)
        mol = geometry.removeAtomFromMol(mol, hidx)
    elif '+' in adduct:
        atom = adduct.split('+')[-1]
        if atom == 'Na':
            mol = geometry.addAtomToMol(mol, atom, idx, covalent=False)
        elif atom == 'H':
            mol = geometry.addAtomToMol(mol, atom, idx, covalent=True)

    _builder = ob.OBBuilder()
    _builder.Build(mol.OBMol)
    mol.localopt(forcefield=forcefield, steps=steps)

    # # adjust partial charge
    # if adduct == '-H':
    #     mol.atoms[idx].OBAtom.SetPartialCharge(mol.atoms[idx].partialcharge - 1)

    return mol


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate adduct from base geometry')
    parser.add_argument('infile', help='path to input .mol2 file')
    parser.add_argument('pkafile', help='path to input .pka file')
    parser.add_argument('adduct', help='adduct type')
    parser.add_argument('mol2', help='path to output .mol2 file')
    parser.add_argument('xyz', help='path to output .xyz file')
    parser.add_argument('pdb', help='path to output .pdb file')
    parser.add_argument('charge', help='path to output .charge file')
    parser.add_argument('-ff', '--forcefield', type=str, default='gaff', help='forcefield type')
    parser.add_argument('-s', '--steps', type=int, default=500, help='number of forcefield optimization steps')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    args = parser.parse_args()
    args.adduct = args.adduct[1:-1]

    # read input
    mol = read_mol(args.infile, fmt='mol2')

    if args.adduct == 'neutral':
        pass
    else:
        pka = read_pka(args.pkafile)

        # generate adduct
        if '+' in args.adduct:
            mol = create_adduct(mol, args.adduct, pka['b1'],
                                forcefield=args.forcefield,
                                steps=args.steps)
        elif '-' in args.adduct:
            mol = create_adduct(mol, args.adduct, pka['a1'],
                                forcefield=args.forcefield,
                                steps=args.steps)

    mol.write('mol2', args.mol2, overwrite=True)
    mol.write('xyz', args.xyz, overwrite=True)
    mol.write('pdb', args.pdb, overwrite=True)

    # reassign charge
    if args.adduct in ['+H', '+Na']:
        # tmp = splitext(args.mol2)[0] + '.tmp.mol2'
        # subprocess.call('obabel %s -O %s --partialcharge eem' % (args.mol2, tmp), shell=True)
        # subprocess.call('mv %s %s' % (tmp, args.mol2), shell=True)
        write_string('1.0', args.charge)

    elif args.adduct == '-H':
        write_string('-1.0', args.charge)

    elif args.adduct == 'neutral':
        write_string('0.0', args.charge)
