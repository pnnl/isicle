import argparse
from string import Template
from pkg_resources import resource_filename
from isicle import __version__


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Prepare molecule for molecular dynamics')
    parser.add_argument('mol2', help='path to .mol2 file')
    parser.add_argument('outfile', help='path to output .mdin file')
    parser.add_argument('-v', '--version', action='version', version=__version__, help='print version and exit')

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--em', action='store_true', help='energy minimization mode')
    mode.add_argument('--iter0', action='store_true', help='first simulated annealing iteration mode')
    mode.add_argument('--sa', action='store_true', help='subsequent simulated annealing iteration mode')

    args = parser.parse_args()

    if args.em is True:
        with open(resource_filename('isicle', 'resources/amber/sander_em.template'), 'r') as f:
            t = Template(f.read())

        d = {'mol2': args.mol2,
             'imin': 1,
             'maxcyc': 500,
             'ncyc': 250,
             'nscm': 1,
             'ntb': 0,
             'igb': 0,
             'cut': 999,
             'rgbmax': 999}

        with open(args.outfile, 'w') as f:
            f.write(t.substitute(d))

    elif args.iter0 is True:
        with open(resource_filename('isicle', 'resources/amber/sander_md0.template'), 'r') as f:
            t = Template(f.read())

        d = {'mol2': args.mol2,
             'imin': 0,
             'ntb': 0,
             'ntf': 2,
             'ntc': 2,
             'ntx': 1,
             'igb': 0,
             'ntpr': 100,
             'ntwx': 100,
             'ntt': 3,
             'nscm': 1,
             'gamma_ln': 1.0,
             'tempi': 0.0,
             'temp0': 300.0,
             'nstlim': 100000,
             'dt': 0.0005,
             'cut': 999,
             'ig': -1}

        with open(args.outfile, 'w') as f:
            f.write(t.substitute(d))

    elif args.sa is True:
        with open(resource_filename('isicle', 'resources/amber/sander_anneal.template'), 'r') as f:
            t = Template(f.read())

        d = {'mol2': args.mol2}

        with open(args.outfile, 'w') as f:
            f.write(t.substitute(d))
