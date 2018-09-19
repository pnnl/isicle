from os.path import *
from string import Template
import argparse


class NWChemHelper:
    def __init__(self, file):
        self.file = file
        self.dir = dirname(self.file)

    def generateNWChemFile(self, template):
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
            elif "+de" in d['filename'].lower():
                d['charge'] = -1
            elif "+ne" in d['filename'].lower():
                d['charge'] = 0

            outfile = splitext(self.file)[0] + '.nw'
            with open(outfile, 'w') as outf:
                outf.write(orig.substitute(d))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Perform the proper preparation steps for an NWChem job.')
    parser.add_argument('file', help='Path to input .xyz file.')
    parser.add_argument('--template', help='Path to template .nw file.', default='-1')

    args = parser.parse_args()

    # if template not specified, pass default template
    if args.template == '-1':
        args.template = 'resources/nwchem/template.nw'

    nwc = NWChemHelper(args.file)
    nwc.generateNWChemFile(args.template)
