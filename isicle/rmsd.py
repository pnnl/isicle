import pybel
import openbabel
import argparse


__version__ = '0.1.0'


def rmsd(mol1, mol2):
    a = next(pybel.readfile("xyz", mol1))
    b = next(pybel.readfile("xyz", mol2))

    align = openbabel.OBAlign(False, True)

    align.SetRefMol(a.OBMol)
    align.SetTargetMol(b.OBMol)
    align.Align()
    return align.GetRMSD()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate RMSD among molecules.')
    parser.add_argument('ref', help='Reference .xyz file.')
    parser.add_argument('outfile', help='Path to output .rmsd file.')
    parser.add_argument('infiles', nargs='+', help='Input .xyz files.')

    args = parser.parse_args()

    total = 0
    for mol in args.infiles:
        total += rmsd(args.ref, mol)

    with open(args.outfile, 'w') as f:
        f.write(str(total) + '\n')
