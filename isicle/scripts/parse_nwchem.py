import pandas as pd
from os.path import *
import glob
from pkg_resources import resource_filename
import argparse


__version__ = '0.1.0'


def parse_dft(path):
    with open(path, 'r') as f:
        lines = f.readlines()

    energy = []
    charges = []
    ready = False
    for line in lines:
        if 'Total DFT energy' in line:
            energy.append(float(line.split()[-1]))

        elif 'Atom       Charge   Shell Charges' in line:
            ready = True
            charges = []
        elif ready is True and line.strip() == '':
            ready = False
        elif ready is True:
            charges.append(line)

    # grab last energy
    energy = energy[-1]

    # process charge information
    df = pd.DataFrame([x.split()[0:4] for x in charges[1:]],
                      columns=['idx', 'Atom', 'Number', 'Charge'])
    df.Number = df.Number.astype('int')
    df.Charge = df.Number - df.Charge.astype('float')

    return energy, df.Charge.values


def extract_geometry(path):
    search = splitext(path)[0]
    geoms = glob.glob(search + '*.xyz')

    if len(geoms) < 1:
        raise IOError("No geometry files found.")

    geoms.sort()

    return geoms[-1]


def generate_mfj(xyz, charges, outfile, masses=resource_filename('isicle', 'resources/mobcal/atomic_mass.tsv')):
    mass = pd.read_csv(masses, sep='\t', usecols=['Number', 'Mass'])

    data = pd.read_csv(xyz, skiprows=2, header=None, delim_whitespace=True, names=['Atom', 'x', 'y', 'z'])
    data['Charge'] = charges

    # merge with atomic masses
    data = pd.merge(data, mass)
    data = data[['x', 'y', 'z', 'Mass', 'Charge']]

    with open(outfile, 'w') as f:
        f.write(outname + '\n')
        f.write('1\n')
        f.write(str(natoms) + '\n')
        f.write('ang\n')
        f.write('calc\n')
        f.write('1.000\n')

        for row in data.values:
            f.write('\t'.join([str(x) for x in row]) + '\n')


def parse_shielding(path, outfile):
    with open(path, 'r') as f:
        lines = f.readlines()

    energy = []
    shield_values = []
    ready = False
    for line in lines:
        if "Total DFT energy" in line:
            energy.append(float(line.split()[-1]))
        elif "Atom:" in line:
            idx = int(line.split()[1])
            atom = line.split()[2]
            ready = True
        elif "isotropic" in line and ready is True:
            shield = float(line.split()[-1])
            shield_values.append([idx, atom, shield])
            ready = False
        elif 'SHIELDING' in line:
            true_idx = [int(x) for x in line.split()[2:]]

    df = pd.DataFrame(shield_values, columns=['index', 'atom', 'shielding'])
    df['dft_energy'] = energy[-1]
    df['index'] = true_idx

    df.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse NWChem DFT output.')
    parser.add_argument('infile', help='Path to NWChem .out file.')
    parser.add_argument('outfile', help='Path to output file.')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='Print version and exit.')

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--dft', action='store_true', help='DFT mode.')
    mode.add_argument('--shielding', action='store_true', help='Shielding mode.')

    args = parser.parse_args()

    if args.dft is True:
        energy, charges = parse_dft(args.infile)
        xyz = extract_geometry(args.infile)
        generate_mfj(xyz, charges, args.outfile)

        # write .energy file
        with open(splitext(outfile)[0] + '.energy', 'w') as f:
            f.write(str(energy))

    elif args.shielding is True:
        parse_shielding(args.infile, args.outfile)
