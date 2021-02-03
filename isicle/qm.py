'''
Assumes NWChem as default quantum chemistry package. Other options can be added later:
 Orca, Quantum Espresso, Gaussian, Schrodinger, etc.
'''

from isicle.interfaces import QMWrapperInterface
from isicle.geometry import load
from isicle.parse import NWChemParser
from isicle.utils import safelist
import tempfile
import os
from string import Template
import itertools
import subprocess
from itertools import cycle


def dft(self, program='NWChem', functional='b3lyp', basis_set='6-31g*'):
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
        task_map = {'optimize': self._configure_optimize,
                    'shielding': self._configure_shielding,
                    'spin': self._configure_spin}

    def load_geometry(self, path):
        # Load geometry
        self.geom = load(path)

        # Extract filename
        self.geom.basename = os.path.splitext(self.geom.path)[0]

    def set_geometry(self, geom):
        # Assign geometry
        self.geom = geom

        # Extract filename
        self.geom.basename = os.path.splitext(self.geom.path)[0]

    def save_geometry(self, fmt='xyz'):
        # Save to temporary directory
        self.fmt = fmt.lower()
        outfile = os.path.join(self.temp_dir, '{}.{}'.format(self.geom.basename, self.fmt.lower()))
        self.geom.save(outfile)

    def _atom_indices(self, atoms=['C', 'H'],
                      lookup={'C': 6, 'H': 1, 'N': 7, 'O': 8, 'F': 9, 'P': 15}):
        atoms = [lookup[x] for x in atoms]
        idx = []
        for a in self.geom.mol.atoms:
            if a.atomicnum in atoms:
                idx.append(a.idx)

        return idx

    def _configure_header(scratch_dir='/scratch', mem_global=1600, mem_heap=100, mem_stack=600):
        d = {'basename': self.geom.basename
             'self.fmt': self.fmt,
             'dirname': self.temp_dir,
             'mem_global': mem_global,
             'mem_heap': mem_heap,
             'mem_stack': mem_stack,
             'scratch_dir': scratch_dir}

        return ('title "{basename}.{fmt}"\n'
                'start {basename}.{fmt}\n\n'
                'memory global {mem_global} mb heap {mem_heap} mb stack {mem_stack} mb\n\n'
                'permanent_dir {dirname}\n'
                'scratch_dir {scratch_dir}\n\n'
                'echo\n'
                'print low\n').format(**d)

    def _configure_load(self, charge=0):
        d = {'basename': self.geom.basename,
             'fmt': self.fmt,
             'dirname': self.temp_dir,
             'charge': charge}

        return ('\ncharge {charge}\n'
                'geometry noautoz noautosym\n'
                ' load {dirname}/{basename}.{fmt}\n'
                'end\n').format(**d)

    def _configure_basis(self, basis_set='6-31G*', ao_basis='cartesian'):

        d = {'ao_basis': ao_basis,
             'basis_set': basis_set}

        return ('\nbasis {ao_basis}\n'
                ' * library {basis_set}\n'
                'end\n').format(**d)

    def _configure_dft(self, functional='b3lyp'):

        d = {'functional': functional}
        
        return ('\ndft\n'
                ' direct\n'
                ' xc {functional}\n'
                ' mulliken\n'               # Do we need this line?
                ' print "mulliken ao"\n'    # (and this one?)
                'end\n').format(**d)

    def _configure_driver(self, max_iter=150):

        d= {'basename': self.geom.basename,
            'fmt': self.fmt,
            'max_iter': max_iter}

        return ('\ndriver\n'
                ' max_iter {max_iter}\n'
                ' {fmt} {basename}_geom\n'
                'end\n')

    def _configure_cosmo(self, solvent='H20', gas=False):

        d = {'solvent': solvent,
             'gas': gas}

        return ('\ncosmo\n'
                ' do_gasphase {gas}\n'
                ' solvent {solvent}\n'
                'end\n').format(**d)

    def _configure_frequency(self, temp=298.15):
        return ('\nfreq'
                ' temp 1 {temp}\n'.format(temp)
                'end\n')

    def _configure_optimize(self, basis_set='6-31G*', ao_basis='cartesian',
                            functional='b3lyp', max_iter=150,
                            cosmo=False, solvent='H20', gas=False,
                            frequency=True, temp=298.15, **kwargs):
        # Add basis block
        s = self._configure_basis(basis_set=basis_set, ao_basis=ao_basis)

        # Add DFT block
        s += self._configure_dft(functional=functional)

        # Add driver block
        s += self._configure_driver(max_iter=150)

        # Add COSMO block
        if cosmo:
            s += self._configure_cosmo(solvent=solvent, gas=gas)

        # Add frequency block
        if frequency:
            s += self._configure_frequency(temp=temp)

        # Add optimize task
        s += '\ntask dft optimize ignore\n'

        # Add frequency task
        if frequency:
            s += 'task dft frequency ignore\n'

        return s

    def _configure_shielding(self, idx, basis_set='6-31G*', ao_basis='cartesian',
                             functional='b3lyp', cosmo=True, solvent='H20', gas=False, energy=True, **kwargs):
        # Add basis block
        s = self._configure_basis(basis_set=basis_set, ao_basis=ao_basis)

        # Add DFT block
        s += self._configure_dft(functional=functional)

        # Add COSMO block
        if cosmo:
            s += self._configure_cosmo(solvent=solvent, gas=gas)

        # Add shielding block
        d = {'ncount': len(idx),
             'nuclei': ' '.join(['%s' % x for x in idx])}

        s += ('\nproperty\n'
              ' SHIELDING {ncount} {nuclei}\n'
              'end\n').format(**d)

        # Add energy task
        if energy:
            s += '\ntask dft energy ignore\n'

        # Add property task
        s += 'task dft property ignore\n'

        return s

    def _configure_spin(self, idx, max_pairs=30, basis_set='6-31G*', ao_basis='cartesian',
                        functional='b3lyp', cosmo=True, solvent='H20', gas=False, energy=True, **kwargs):
        # Enumerate spin-spin couplings
        pairs = list(itertools.combinations(idx, 2))

        # Add basis block
        s = self._configure_basis(basis_set=basis_set, ao_basis=ao_basis)

        # Add DFT block
        s += self._configure_dft(functional=functional)

        # Add COSMO block
        if cosmo:
            s += self._configure_cosmo(solvent=solvent, gas=gas)

        # Add spin block
        for i, p in enumerate(pairs):
            if i % max_pairs == 0:
                if i > 0:
                    s += '\nend\n'
                    s += '\ntask dft property\n\n'
                s += 'property\n'
                npairs = min(len(pairs) - i, max_pairs)
                s += ' SPINSPIN {}}'.format(npairs)
            s += ' {} {}'.format(*p)

        s += '\nend\n'

        # Add energy task
        if energy:
            s += '\ntask dft energy ignore\n'

        # Add property task
        s += 'task dft property ignore\n'

        return s

    def configure(self, tasks='optimize', functional='b3lyp', basis_set='6-31g*', ao_basis='cartesian',
                  charge=0, atoms=['C', 'H'], energy=True, frequency=True, temp=298.15,
                  cosmo=False, solvent='H20', gas=False, max_iter=150,
                  mem_global=1600, mem_heap=100, mem_stack=600, scratch_dir='/scratch'):

        # Cast to list safely
        tasks = safelist(tasks)
        functional = safelist(functional)
        basis_set = safelist(basis_set)
        ao_basis = safelist(ao_basis)
        max_iter = safelist(max_iter)
        atoms = safelist(atoms)
        cosmo = safelist(cosmo)

        # Container for final configuration script
        nwconfig = ''

        # Check lengths
        if not ((len(tasks) == len(functional)) or (len(functional) == 1)):
            raise ValueError('Functional must be assigned globally or per task.')

        if not ((len(tasks) == len(basis_set)) or (len(basis_set) == 1)):
            raise ValueError('Basis set must be assigned globally or per task.')

        if not ((len(tasks) == len(ao_basis)) or (len(ao_basis) == 1)):
            raise ValueError('AO basis must be assigned globally or per task.')

        if not ((len(tasks) == len(cosmo)) or (len(cosmo) == 1)):
            raise ValueError('Maximum iterations must be assigned globally or per task.')

        # Extract atom index information
        idx = self._atom_indices(atoms=atoms)

        # Generate header information
        nwconfig += self._configure_header(scratch_dir=scratch_dir, mem_global=mem_global, mem_heap=mem_heap, mem_stack=mem_stack)

        # Load geometry
        nwconfig += self._configure_load(charge=charge)

        # Configure tasks
        for task, f, b, a, c in zip(tasks, cycle(functional), cycle(basis), cycle(ao_basis), cycle(cosmo)):
            # TODO: finish this
            nwconfig += self.task_map[task](functional=f, basis_set=b, ao_basis=a,
                                            energy=energy, frequency=frequency, temp=temp,
                                            cosmo=c, solvent=solvent, gas=gas, max_iter=max_iter)

        # Write to file
        with open(os.path.join(self.temp_dir, self.geom.basename + '.nw'), 'w') as f:
            f.write(nwconfig)

    def run(self):
        infile = os.path.join(self.temp_dir, self.geom.basename + '.nw')
        outfile = os.path.join(self.temp_dir, self.geom.basename + '.out')
        logfile = os.path.join(self.temp_dir, self.geom.basename + '.log')
        subprocess.call('nwchem {} > {} 2> {}'.format(infile, outfile, logfile, shell=True))

    def finish(self, keep_files=True, path=None):
        parser = NWChemParser()
        parser.load(os.path.join(self.temp_dir, self.geom.basename + '.out'))
        result = parser.parse(to_parse=self.to_parse)

        if keep_files is True:
            import shutil
            import glob

            if path is None:
                raise ValueError('Must supply `path`.') 
            else:
                # TODO: anything else to keep?
                shutil.copy2(os.path.join(self.temp_dir, self.geom.basename + '.nw'), path)
                shutil.copy2(os.path.join(self.temp_dir, self.geom.basename + '.out'), path)
                shutil.copy2(os.path.join(self.temp_dir, self.geom.basename + '.log'), path)

                geoms = glob.glob(self.temp_dir, '*.{}'.format(self.fmt))
                [shutil.copy2(x, path) for x in geoms]

        # Remove temporary files
        os.removedirs(self.temp_dir)

        return result

# TODO: add other DFT methods as needed
