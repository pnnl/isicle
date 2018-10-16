import argparse


__version__ = '0.1.0'


class NWChemHelper:
    def __init__(self, file):
        self.file = file
        self.dir = dirname(self.file)

    def dft(self, template):
        with open(template, 'r') as t:
            orig = Template(t.read())

        d = {'filename': None, 'dir': self.dir, 'charge': None}
        d['filename'] = splitext(basename(self.file))[0]

        if not isfile(d['filename']):
            if "+2h" in d['filename'].lower():
                d['charge'] = 2
            elif "+3h" in d['filename'].lower():
                d['charge'] = 3
            elif "+2na" in d['filename'].lower():
                d['charge'] = 2
            elif "+h+na" in d['filename'].lower():
                d['charge'] = 2
            elif "+h" in d['filename'].lower():
                d['charge'] = 1
            elif "+na" in d['filename'].lower():
                d['charge'] = 1
            elif "-h" in d['filename'].lower():
                d['charge'] = -1
            elif "_ne_" in d['filename'].lower():
                d['charge'] = 0

            outfile = splitext(self.file)[0] + '.nw'
            with open(outfile, 'w') as outf:
                outf.write(orig.substitute(d))

    def shielding(self, template, atoms):
        with open(template, 'r') as t:
            orig = Template(t.read())

        d = {'filename': None, 'dir': self.dir, 'charge': None}
        d['filename'] = splitext(basename(self.file))[0]

        idx = self.nuclei(atoms)
        d['ncount'] = len(idx)
        d['nuclei'] = ' '.join(['%s' % x for x in idx])

        outfile = splitext(self.file)[0] + '.nw'
        with open(outfile, 'w') as outf:
            outf.write(orig.substitute(d))

    def nuclei(self, atoms):
        mol = next(pybel.readfile("xyz", args.file))
        idx = []
        for a in mol.atoms:
            if a.atomicnum in atoms:
                idx.append(a.idx)
        return idx


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform the proper preparation steps for an NWChem job.')
    parser.add_argument('file', help='Path to input .xyz file.')
    parser.add_argument('--template', help='Path to template .nw file.', default='default')
    parser.add_argument('--version', '-v', action='version', version=__version__, help='Print version and exit.')

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--dft', action='store_true', help='DFT mode.')
    mode.add_argument('--shielding', action='store_true', help='Shielding mode.')

    parser.add_argument('--shifts', '-s', nargs='+', help='Atomic numbers of atoms to perform shielding calcs.')

    args = parser.parse_args()

    from os.path import *
    from string import Template
    import pybel
    from pkg_resources import resource_filename

    nwc = NWChemHelper(args.file)

    if args.dft is True:
        if args.template == 'default':
            args.template = resource_filename('isicle', 'resources/nwchem/dft.template')

        nwc.dft(args.template)

    elif args.shielding is True:
        if args.shifts is not None:
            lookup = {'C': 6, 'H': 1}
            shifts = [lookup[x] for x in args.shifts]
        else:
            parser.error('Please select atoms to perform shielding calcs.')

        if args.template == 'default':
            args.template = resource_filename('isicle', 'resources/nwchem/shielding.template')

        nwc.shielding(args.template, shifts)
