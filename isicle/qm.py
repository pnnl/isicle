'''
Assumes NWChem as default quantum chemistry package. Other options can be added later:
 Orca, Quantum Espresso, Gaussian, Schrodinger, etc.
'''

from isicle.interfaces import QMWrapperInterface
from isicle.geometry import load
from isicle.parse import NWChemParser
import tempfile
import os
from string import Template
import itertools
import subprocess


def dft(self, program='NWChem', functional='b3lyp', basisset='6-31g*'):
    '''
    Optimize geometry, either XYZ or PDB, using stated functional and basis set.
    Additional inputs can be grid size, optimization criteria level,
    '''
    # TODO: define input, arguments, and output
    # Should return instance (or list) of DFT

    # save geometry to input file
    # load/generate .nw script
    # submit job
    # parse result
    # return decorated `*OptimizedGeometry` class instance
    raise NotImplementedError


class NWChemWrapper(QMWrapperInterface):

    def __init__(self):
        self.temp_dir = tempfile.TemporaryDirectory()
    
    # TODO: what if user wants to pass a `Geometry` instance instead?
    def load_geometry(self, path):
        # Load geometry
        self.geom = load(path)

        # Extract filename
        self.basename = os.path.splitext(os.path.basename(path))[0]

        # Save to temporary directory
        outfile = os.path.join(self.temp_dir, self.basename + '.xyz')
        self.geom.save(outfile)

    def _set_template(self, mode='optimize', path=None):
        # Possible selections
        options = ['optimze', 'shielding', 'spin', 'custom']
        
        # Mode switch
        if mode.lower() == 'optimize':
            path = 'path/to/optimize.template'
            self.to_parse = ['geometry', 'energy']

        elif mode.lower() == 'shielding':
            path = 'path/to/shielding.template'
            self.to_parse = ['geometry', 'shielding', 'energy']

        elif mode.lower() == 'spin':
            path = 'path/to/spin.template'
            self.to_parse = ['geometry', 'shielding', 'energy']

        elif mode.lower() == 'custom':
            if path is None:
                raise ValueError('Must supply `path` for custom template.')

        else:
            raise ValueError('Mode "{}" not recognized. Please supply one of {}.'.format(mode, options))

        # Load template
        with open(path, 'r') as template:
            return Template(template.read())

    def configure(self, mode='optimize', path=None, functional='b3lyp', basisset='6-31g*',
                  atoms=['C', 'H'], solvent='H20'):
        # Load template
        template = self._set_template(mode=mode, path=path)

        # Extract atom index information
        lookup = {'C': 6, 'H': 1}
        atoms = [lookup[x] for x in atoms]
        idx = []
        for a in self.geom.mol.atoms:
            if a.atomicnum in atoms:
                idx.append(a.idx)

        # Enumerate spin-spin couplings
        pairs = list(itertools.combinations(idx, 2))
        s = ''
        maxpairs = 30
        for i, p in enumerate(pairs):
            if i % maxpairs == 0:
                if i > 0:
                    s += '\nend\n'
                    s += '\ntask dft property\n\n'
                s += 'property\n'
                npairs = min(len(pairs) - i, maxpairs)
                s += ' SPINSPIN %i' % npairs
            s += ' %i %i' % p
        s += '\nend\n'
        s += '\ntask dft property'

        # Substitions
        # TODO: sanity check on selections?
        d = {'filename': self.basename + '.xyz',
             'dir': self.temp_dir,
             'functional': functional.lower(),
             'basisset': basisset.lower(),
             'ncount': len(idx),
             'nuclei': ' '.join(['%s' % x for x in idx]),
             'solvent': solvent.lower(),
             'spin': s}

        # Save to temporary directory
        outfile = os.path.join(self.temp_dir, self.basename + '.nw')
        with open(outfile, 'w') as f:
            f.write(orig.substitute(d))

    def run(self):
        infile = os.path.join(self.temp_dir, self.basename + '.nw')
        outfile = os.path.join(self.temp_dir, self.basename + '.out')
        logfile = os.path.join(self.temp_dir, self.basename + '.log')
        subprocess.call('nwchem {} > {} 2> {}'.format(infile, outfile, logfile, shell=True))

    def finish(self, keep_files=True, path=None):
        parser = NWChemParser()
        parser.load(os.path.join(self.temp_dir, self.basename + '.out'))
        result = parser.parse(to_parse=self.to_parse)

        if keep_files is True:
            import shutil
            import glob

            if path is None:
                raise ValueError('Must supply `path`.') 
            else:
                # TODO: anything else to keep?
                shutil.copy2(os.path.join(self.temp_dir, self.basename + '.nw'), path)
                shutil.copy2(os.path.join(self.temp_dir, self.basename + '.out'), path)
                shutil.copy2(os.path.join(self.temp_dir, self.basename + '.log'), path)

                xyzs = glob.glob(self.temp_dir, '*.xyz')
                [shutil.copy2(x, path) for x in xyzs]

        # Remove temporary files
        os.removedirs(self.temp_dir)

        return result

# TODO: add other DFT methods as needed
