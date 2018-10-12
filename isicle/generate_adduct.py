import argparse


__version__ = '0.1.0'


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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate adduct from base geometry.')
    parser.add_argument('infile', help='Path to input .mol file.')
    parser.add_argument('pkafile', help='Path to input .pka file.')
    parser.add_argument('adduct', help='Adduct type.')
    parser.add_argument('mol2', help='Path to output .mol2 file.')
    parser.add_argument('xyz', help='Path to output .xyz file.')
    parser.add_argument('--forcefield', '-ff', type=str, default='gaff', help='Forcefield type.')
    parser.add_argument('--steps', '-s', type=int, default=500, help='Number of forcefield optimization steps.')
    parser.add_argument('--version', '-v', action='version', version=__version__, help='Print version and exit.')

    args = parser.parse_args()

    import pybel
    from isicle.core import geometry
    from isicle.core.utils import read_pka

    args.adduct = args.adduct[1:-1]

    # read inputs
    mol = next(pybel.readfile("mol", args.infile))
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
    else:
        adduct = mol

    adduct.write('mol2', args.mol2, overwrite=True)
    adduct.write('xyz', args.xyz, overwrite=True)
