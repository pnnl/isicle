import argparse
from openbabel import pybel
from isicle import __version__


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
    parser = argparse.ArgumentParser(description='Standardize mol2 for use with OpenBabel')
    parser.add_argument('mol2', help='path to .mol2 file')
    parser.add_argument('ref', help='path to reference .mol2 file')
    parser.add_argument('outfile', help='path to output .xyz file')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    args = parser.parse_args()

    standardize(args.mol2, args.ref, args.outfile)
