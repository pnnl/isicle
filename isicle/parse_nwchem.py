from os.path import *
import glob
import pandas as pd
import shutil
import argparse


__version__ = '0.1.0'


def XYZtoMFJ(resfile, outpath):
    # atomic masses
    masses = pd.read_csv('isicle/resources/mobcal/atomic_mass.tsv', sep='\t', usecols=['Number', 'Mass'])

    # read NWChem output file
    with open(resfile, 'r') as f:
        res = f.readlines()

    # parse output
    natoms, lowdinIdx, energies = parseOutput(res)

    # grab relevant files
    files = glob.glob(splitext(resfile)[0] + '_geom-*.xyz')
    if len(files) < 1:
        raise IOError("No geometry files found.")
    elif len(files) == 1:
        idx = [0]
        geoms = [files[0].rsplit('_', 1)[0] + '_charge.xyz']
        shutil.copyfile(files[0], geom[0])
    else:
        files.sort()

        idx = [0, -1]
        geoms = [files[0].rsplit('_', 1)[0] + '_charge.xyz', files[-1].rsplit('_', 1)[0] + '_geom+charge.xyz']
        shutil.copyfile(files[0], geoms[0])
        shutil.copyfile(files[-1], geoms[1])

    # iterate through successful geometries
    for i, geom in zip(idx, geoms):
        # read xyz file
        xyz = pd.read_csv(geom, skiprows=2, header=None, delim_whitespace=True, names=['Atom', 'x', 'y', 'z'])
        xyz.drop('Atom', 1, inplace=True)

        df = pd.DataFrame([x.split()[0:4] for x in res[lowdinIdx[i]:lowdinIdx[i] + natoms]],
                          columns=['idx', 'Atom', 'Number', 'Charge'])

        df.Number = df.Number.astype('int')
        df.Charge = df.Number - df.Charge.astype('float')

        # merge xyz and NWChem output
        data = pd.merge(xyz, df, left_index=True, right_index=True)

        # merge with atomic masses
        data = pd.merge(data, masses)
        data = data[['x', 'y', 'z', 'Mass', 'Charge']]

        # write mfj file
        outname = basename(splitext(geom)[0] + '.mfj')
        with open(join(outpath, outname), 'w') as f:
            f.write(outname + '\n')
            f.write('1\n')
            f.write(str(natoms) + '\n')
            f.write('ang\n')
            f.write('calc\n')
            f.write('1.000\n')

            for row in data.values:
                f.write('\t'.join([str(x) for x in row]) + '\n')

        # write energy
        ename = splitext(outname)[0] + '.energy'
        with open(join(outpath, ename), 'w') as f:
            f.write(str(energies[i]))


def parseOutput(res, idx=0):
    indices = []
    lowdinIdx = []
    energies = []

    ready = False
    for i, row in enumerate(res):
        if 'No.' in row and len(indices) == 0:
            indices.append(i + 2)  # 0
        elif 'Atomic Mass' in row and len(indices) == 1:
            indices.append(i - 1)  # 1
            indices.append(i + 3)  # 2
        elif 'Effective nuclear repulsion energy' in row and len(indices) == 3:
            indices.append(i - 2)  # 3

        elif 'Lowdin Population Analysis' in row:
            ready = True

        elif 'Shell Charges' in row and ready is True:
            lowdinIdx.append(i + 2)
            ready = False

        if "Total DFT energy" in row:
            energies.append(float(row.rstrip().split('=')[-1]))

    natoms = int(res[indices[1] - 1].split()[0])
    return natoms, lowdinIdx, energies


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse NWChem DFT output.')
    parser.add_argument('infile', help='Path to NWChem .out file.')
    parser.add_argument('outdir', help='Path to output directory.')
    parser.add_argument('--mode', '-m', default='dft', help='Specify NWChem output modality [dft, shielding].')

    args = parser.parse_args()

    if args.mode.lower() == 'dft':
        XYZtoMFJ(args.infile, args.outdir)
