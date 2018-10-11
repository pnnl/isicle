import argparse


__version__ = '0.1.0'


def standardize(path, reference, output):
    traj = next(pybel.readfile("mol2", path))
    ref = next(pybel.readfile("mol2", reference))
    for iatom in traj.atoms:
        ob = iatom.OBAtom
        idx = ob.GetIndex()

        jatom = ref.OBMol.GetAtomById(idx)
        ob.SetType(jatom.GetType())
        ob.SetAtomicNum(jatom.GetAtomicNum())

    traj.write("xyz", output, True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Standardize mol2 for use with OpenBabel.')
    parser.add_argument('mol2', help='Path to .mol2 file.')
    parser.add_argument('ref', help='Path to reference .mol2 file.')
    parser.add_argument('outfile', help='Path to output .xyz file.')
    parser.add_argument('--version', action='version', version=__version__, help='Print version and exit.')

    args = parser.parse_args()

    import pybel

    standardize(args.mol2, args.ref, args.outfile)
