import argparse


__version__ = '0.1.0'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate tleap config.')
    parser.add_argument('mol2', help='Path to .mol2 file.')
    parser.add_argument('frcmod', help='Path to .frcmod file.')
    parser.add_argument('outfile', help='Path to output .config file.')
    parser.add_argument('--version', '-v', action='version', version=__version__, help='Print version and exit.')

    args = parser.parse_args()

    from string import Template
    from os.path import *
    from pkg_resources import resource_filename

    with open(resource_filename('isicle', 'resources/amber/tleap.template'), 'r') as f:
        t = Template(f.read())

    d = {'gaff': resource_filename('isicle', 'resources/amber/leaprc.gaff'),
         'ff': resource_filename('isicle', 'resources/amber/leaprc.ff14SB'),
         'frcmod': args.frcmod,
         'mol2': args.mol2,
         'prmtop': splitext(args.outfile)[0] + '.top',
         'inpcrd': splitext(args.outfile)[0] + '.crd',
         'log': splitext(args.outfile)[0] + '.main.log'}

    with open(args.outfile, 'w') as f:
        f.write(t.substitute(d))
