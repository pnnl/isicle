import argparse
import pybel
import openbabel as ob
from isicle.resources import geometry
from isicle.utils import read_pka
from isicle import __version__


def create_adduct(mol, adduct, idx, forcefield='mmff94', steps=500):
    if adduct == '-H':
        hidx = geometry.nearestHydrogen(mol, idx)
        adduct = geometry.removeAtomFromMol(mol, hidx)
    elif '+' in adduct:
        atom = adduct.split('+')[-1]
        if atom == 'Na':
            adduct = geometry.addAtomToMol(mol, atom, idx, covalent=False)
        elif atom == 'H':
            adduct = geometry.addAtomToMol(mol, atom, idx, covalent=True)

    _builder = ob.OBBuilder()
    _builder.Build(adduct.OBMol)
    adduct.localopt(forcefield=forcefield, steps=steps)

    # adjust partial charge
    if adduct == '-H':
        adduct.atoms[pka['a1']].OBAtom.SetPartialCharge(adduct.atoms[idx].partialcharge - 1)

    return adduct


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate adduct from base geometry')
    parser.add_argument('infile', help='path to input .mol file')
    parser.add_argument('pkafile', help='path to input .pka file')
    parser.add_argument('adduct', help='adduct type')
    parser.add_argument('mol2', help='path to output .mol2 file')
    parser.add_argument('xyz', help='path to output .xyz file')
    parser.add_argument('-ff', '--forcefield', type=str, default='gaff', help='forcefield type')
    parser.add_argument('-s', '--steps', type=int, default=500, help='number of forcefield optimization steps')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    args = parser.parse_args()
    args.adduct = args.adduct[1:-1]

    # read inputs
    mol = next(pybel.readfile("mol", args.infile))

    if args.adduct == 'neutral':
        adduct = mol
    else:
        pka = read_pka(args.pkafile)

        # generate adduct
        if '+' in args.adduct:
            adduct = create_adduct(mol, args.adduct, pka['b1'],
                                   forcefield=args.forcefield,
                                   steps=args.steps)
        elif '-' in args.adduct:
            adduct = create_adduct(mol, args.adduct, pka['a1'],
                                   forcefield=args.forcefield,
                                   steps=args.steps)

    adduct.write('mol2', args.mol2, overwrite=True)
    adduct.write('xyz', args.xyz, overwrite=True)
