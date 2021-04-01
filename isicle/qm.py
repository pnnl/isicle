from isicle.interfaces import WrapperInterface
from isicle.parse import NWChemParser
from isicle.utils import safelist
import tempfile
import os
from string import Template
from itertools import combinations, cycle
import subprocess


def _program_selector(program):
    '''
    Selects a supported quantum mechanical program for associated simulation.
    Currently only NWChem has been implemented.

    Parameters
    ----------
    program : str
        Alias for program selection (e.g. NWChem).

    Returns
    -------
    program
        Wrapped functionality of the selected program. Must implement
        :class:`~isicle.interfaces.QMWrapperInterface`.

    '''

    program_map = {'nwchem': NWChemWrapper}

    if program.lower() in program_map.keys():
        return program_map[program.lower()]()
    else:
        raise ValueError(('{} not a supported quantum mechanical program.')
                         .format(program))


def dft(geom, program='NWChem', template=None, **kwargs):
    '''
    Optimize geometry via density functional theory using supplied functional
    and basis set.

    Parameters
    ----------
    geom : :obj:`~isicle.geometry.Geometry`
        Molecule representation.
    program : str
        Alias for program selection (NWChem).
    template : str
        Path to optional template to bypass default configuration process.
    **kwargs
        Keyword arguments to configure the simulation.
        See :meth:`~isicle.qm.NWChemWrapper.configure`.

    Returns
    -------
    :obj:`~isicle.geometry.Geometry`
        DFT-optimized molecule representation.
    :obj:`~isicle.parse.NWChemResult`
        Result object containing relevant outputs from the simulation.

    '''

    # Select program
    qmw = _program_selector(program)

    # Set geometry
    qmw.set_geometry(geom)

    # Save geometry
    qmw.save_geometry(fmt=kwargs.pop('fmt'))

    # Configure
    if template is not None:
        qmw.configure_from_template(template)
    else:
        qmw.configure(**kwargs)

    # Save configuration file
    qmw.save_config()

    # Run QM simulation
    qmw.run()

    # Finish/clean up
    res = qmw.finish()

    # Create new Geometry with updated structure
    # res['geometry'] will be None or a path to an xyz file.
    geom = geom._update_structure(False, xyz_filename=res['geometry'])

    # Erase old properties and add new event and DFT properties
    geom.global_properties = {}
    geom._update_history('dft')
    geom = geom.add_global_properties(res.to_dict())

    return geom, res


class NWChemWrapper(WrapperInterface):
    '''
    Wrapper for NWChem functionality.

    Implements :class:`~isicle.interfaces.QMWrapperInterface` to ensure
    required methods are exposed.

    Attributes
    ----------
    temp_dir : str
        Path to temporary directory used for simulation.
    task_map : dict
        Alias mapper for supported quantum mechanical presets. Thses include
        "optimze", "shielding", and "spin".
    geom : :obj:`~isicle.geometry.Geometry`
        Internal molecule representation.
    fmt : str
        File extension indicator.
    config : str
        Configuration information for simulation.

    '''

    def __init__(self):
        '''
        Initialize :obj:`~isicle.qm.NWChemWrapper` instance.

        Creates temporary directory for intermediate files, establishes aliases
        for preconfigured tasks.

        '''

        self.temp_dir = tempfile.TemporaryDirectory()
        self.task_map = {'optimize': self._configure_optimize,
                         'shielding': self._configure_shielding,
                         'spin': self._configure_spin}

    def set_geometry(self, geom):
        '''
        Set :obj:`~isicle.geometry.Geometry` instance for simulation.

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.

        '''

        # Assign geometry
        self.geom = geom

    def save_geometry(self, fmt='xyz'):
        '''
        Save internal :obj:`~isicle.geometry.Geometry` representation to file.

        Parameters
        ----------
        fmt : str
            Filetype used by NWChem. Must be "xyz" or "pdb."

        Raises
        ------
        TypeError
            If geometry loaded from .xyz is saved to another format.

        '''

        # Path operations
        self.fmt = fmt.lower()
        outfile = os.path.join(self.temp_dir.name,
                               '{}.{}'.format(self.geom.basename,
                                              self.fmt.lower()))

        # All other formats
        self.geom.save(outfile)

    def _configure_header(self, scratch_dir='/scratch', mem_global=1600,
                          mem_heap=100, mem_stack=600):
        '''
        Generate header block of NWChem configuration.

        Parameters
        ----------
        scratch_dir : str
            Path to simulation scratch directory.
        mem_global : int
            Global memory allocation in MB.
        mem_heap : int
            Heap memory allocation in MB.
        mem_stack : int
            Stack memory allocation in MB.

        Returns
        -------
        str
            Header block of NWChem configuration.

        '''

        d = {'basename': self.geom.basename,
             'dirname': self.temp_dir.name,
             'mem_global': mem_global,
             'mem_heap': mem_heap,
             'mem_stack': mem_stack,
             'scratch_dir': scratch_dir}

        return ('title "{basename}"\n'
                'start {basename}\n\n'
                'memory global {mem_global} mb heap {mem_heap} '
                'mb stack {mem_stack} mb\n\n'
                'permanent_dir {dirname}\n'
                'scratch_dir {scratch_dir}\n\n'
                'echo\n'
                'print low\n').format(**d)

    def _configure_load(self, charge=0):
        '''
        Generate geometry load block of NWChem configuration.

        Parameters
        ----------
        charge : int
            Nominal charge of the molecule to be optimized.

        Returns
        -------
        str
            Geometry load block of NWChem configuration.

        '''

        d = {'basename': self.geom.basename,
             'fmt': self.fmt,
             'dirname': self.temp_dir.name,
             'charge': charge}

        return ('\ncharge {charge}\n'
                'geometry noautoz noautosym\n'
                ' load {dirname}/{basename}.{fmt}\n'
                'end\n').format(**d)

    def _configure_basis(self, basis_set='6-31G*', ao_basis='cartesian'):
        '''
        Generate basis set block of NWChem configuration.

        Parameters
        ----------
        basis_set : str
            Basis set selection.
        ao_basis : str
            Angular function selection ("spherical", "cartesian").

        Returns
        -------
        str
            Basis set block of NWChem configuration.

        '''

        d = {'ao_basis': ao_basis,
             'basis_set': basis_set}

        return ('\nbasis {ao_basis}\n'
                ' * library {basis_set}\n'
                'end\n').format(**d)

    def _configure_dft(self, functional='b3lyp', odft=False):
        '''
        Generate DFT block of NWChem configuration.

        Parameters
        ----------
        functional : str
            Functional selection.
        odft : bool
            Indicate whether to use open DFT functional (required for spin-spin
            couplings).

        Returns
        -------
        str
            DFT block of NWChem configuration.

        '''

        d = {'functional': functional,
             'dft': 'odft' if odft is True else 'dft'}

        return ('\n{dft}\n'
                ' direct\n'
                ' xc {functional}\n'
                ' mulliken\n'               # Do we need this line?
                ' print "mulliken ao"\n'    # (and this one?)
                'end\n').format(**d)

    def _configure_driver(self, max_iter=150):
        '''
        Generate driver block of NWChem configuration.

        Parameters
        ----------
        max_iter : int
            Maximum number of optimization iterations.

        Returns
        -------
        str
            Driver block of NWChem configuration.

        '''

        d = {'basename': self.geom.basename,
             'fmt': self.fmt,
             'max_iter': max_iter}

        return ('\ndriver\n'
                ' max_iter {max_iter}\n'
                ' {fmt} {basename}_geom\n'
                'end\n').format(**d)

    def _configure_cosmo(self, solvent='H2O', gas=False):
        '''
        Generate COSMO block of NWChem configuration.

        Parameters
        ----------
        solvent : str
            Solvent selection.
        gas : bool
            Indicate whether to use gas phase calculations.

        Returns
        -------
        str
            COSMO block of NWChem configuration.

        '''

        d = {'solvent': solvent,
             'gas': gas}

        return ('\ncosmo\n'
                ' do_gasphase {gas}\n'
                ' solvent {solvent}\n'
                'end\n').format(**d)

    def _configure_frequency(self, temp=298.15):
        '''
        Configure frequency block of NWChem configuration.

        Parameters
        ----------
        temp : float
            Temperature for frequency calculation.

        Returns
        -------
        str
            Frequency block of NWChem configuration.

        '''

        return ('\nfreq\n'
                ' temp 1 {}\n'
                'end\n').format(temp)

    def _configure_optimize(self, basis_set='6-31G*', ao_basis='cartesian',
                            functional='b3lyp', max_iter=150,
                            cosmo=False, solvent='H2O', gas=False,
                            frequency=True, temp=298.15, **kwargs):
        '''
        Generate meta optimization block of NWChem configuration.

        Includes basis, DFT, and driver blocks; can include COSMO and/or
        frequency blocks.

        Parameters
        ----------
        basis_set : str
            Basis set selection.
        ao_basis : str
            Angular function selection ("spherical", "cartesian").
        functional : str
            Functional selection.
        max_iter : int
            Maximum number of optimization iterations.
        cosmo : bool
            Indicate whether to include COSMO block.
        solvent : str
            Solvent selection. Only used if `cosmo` is True.
        gas : bool
            Indicate whether to use gas phase calculations. Only used if
            `cosmo` is True.
        frequency : bool
            Indicate whether to include frequency block.
        temp : float
            Temperature for frequency calculation. Only used if `frequency` is
            True.
        **kwargs
            Arbitrary additional arguments (unused).

        Returns
        -------
        str
            Optimization meta block of NWChem configuration.

        '''

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

    def _configure_shielding(self, basis_set='6-31G*', ao_basis='cartesian',
                             functional='b3lyp', cosmo=True, solvent='H2O',
                             gas=False, energy=True, **kwargs):
        '''
        Generate meta shielding block of NWChem configuration.

        Includes basis and DFT blocks; can include COSMO and/or single-point
        energy calculation blocks.

        Parameters
        ----------
        basis_set : str
            Basis set selection.
        ao_basis : str
            Angular function selection ("spherical", "cartesian").
        functional : str
            Functional selection.
        max_iter : int
            Maximum number of optimization iterations.
        cosmo : bool
            Indicate whether to include COSMO block.
        solvent : str
            Solvent selection. Only used if `cosmo` is True.
        gas : bool
            Indicate whether to use gas phase calculations. Only used if
            `cosmo` is True.
        energy : bool
            Indicate whether to include single-point energy calculation block.
        **kwargs
            Arbitrary additional arguments (unused).

        Returns
        -------
        str
            Shielding meta block of NWChem configuration.

        '''

        # Add basis block
        s = self._configure_basis(basis_set=basis_set, ao_basis=ao_basis)

        # Add DFT block
        s += self._configure_dft(functional=functional)

        # Add COSMO block
        if cosmo:
            s += self._configure_cosmo(solvent=solvent, gas=gas)

        # Add shielding block
        d = {'ncount': len(self.idx),
             'nuclei': ' '.join(['%s' % x for x in self.idx])}

        s += ('\nproperty\n'
              ' SHIELDING {ncount} {nuclei}\n'
              'end\n').format(**d)

        # Add energy task
        if energy:
            s += '\ntask dft energy ignore\n'

        # Add property task
        s += 'task dft property ignore\n'

        return s

    def _configure_spin(self, max_pairs=30, basis_set='6-31G*',
                        ao_basis='cartesian', functional='b3lyp', cosmo=True,
                        solvent='H2O', gas=False, energy=True, **kwargs):
        '''
        Generate meta spin-spin coupling block of NWChem configuration.

        Includes basis and DFT blocks; can include COSMO and/or single-point
        energy calculation blocks.

        Parameters
        ----------
        max_pairs : int
            Maximum number of spin-spin pairs per spin-spin coupling block.
            Note: do not modify.
        basis_set : str
            Basis set selection.
        ao_basis : str
            Angular function selection ("spherical", "cartesian").
        functional : str
            Functional selection.
        cosmo : bool
            Indicate whether to include COSMO block.
        solvent : str
            Solvent selection. Only used if `cosmo` is True.
        gas : bool
            Indicate whether to use gas phase calculations. Only used if
            `cosmo` is True.
        energy : bool
            Indicate whether to include single-point energy calculation block.
        **kwargs
            Arbitrary additional arguments (unused).

        Returns
        -------
        str
            Spin-spin coupling meta block of NWChem configuration.

        '''

        # Enumerate spin-spin couplings
        pairs = list(combinations(self.idx, 2))

        # Add basis block
        s = self._configure_basis(basis_set=basis_set, ao_basis=ao_basis)

        # Add DFT block
        s += self._configure_dft(functional=functional, odft=True)

        # Add COSMO block
        if cosmo:
            s += self._configure_cosmo(solvent=solvent, gas=gas)

        # Add spin block
        for i, p in enumerate(pairs):
            if i % max_pairs == 0:
                if i > 0:
                    s += '\nend\n'
                    s += '\ntask dft property\n\n'
                s += '\nproperty\n'
                npairs = min(len(pairs) - i, max_pairs)
                s += ' SPINSPIN {}'.format(npairs)
            s += ' {} {}'.format(*p)

        s += '\nend\n'

        # Add energy task
        if energy:
            s += '\ntask dft energy ignore\n'

        # Add property task
        s += 'task dft property ignore\n'

        return s

    def configure(self, tasks='optimize', functional='b3lyp',
                  basis_set='6-31g*', ao_basis='cartesian', charge=0,
                  atoms=['C', 'H'], energy=True, frequency=True, temp=298.15,
                  cosmo=False, solvent='H2O', gas=False, max_iter=150,
                  mem_global=1600, mem_heap=100, mem_stack=600,
                  scratch_dir='/scratch'):
        '''
        Configure NWChem simulation.

        Parameters
        ----------
        tasks : str or list of str
            Tasks text.
        functional : str or list of str
            Functional selection. Supply globally or per task.
        basis_set : str or list of str
            Basis set selection. Supply globally or per task.
        ao_basis : str or list of str
            Angular function selection ("spherical", "cartesian"). Supply
            globally or per task.
        charge : int
            Nominal charge of the molecule to be optimized.
        atoms : list of str
            Atom types of interest.
        energy : bool
            Indicate whether to include single-point energy calculation block.
        frequency : bool
            Indicate whether to include frequency block.
        temp : float
            Temperature for frequency calculation. Only used if `frequency` is
            True.
        cosmo : bool
            Indicate whether to include COSMO block. Supply globally or per
            task.
        solvent : str
            Solvent selection. Only used if `cosmo` is True. Supply globally or
            per task.
        gas : bool
            Indicate whether to use gas phase calculations. Only used if
            `cosmo` is True. Supply globally or per task.
        max_iter : int
            Maximum number of optimization iterations.
        scratch_dir : str
            Path to simulation scratch directory.
        mem_global : int
            Global memory allocation in MB.
        mem_heap : int
            Heap memory allocation in MB.
        mem_stack : int
            Stack memory allocation in MB.

        Returns
        -------
        str
            NWChem configuration.

        '''

        # Cast to list safely
        tasks = safelist(tasks)
        functional = safelist(functional)
        basis_set = safelist(basis_set)
        ao_basis = safelist(ao_basis)
        atoms = safelist(atoms)
        cosmo = safelist(cosmo)

        # Container for final configuration script
        config = ''

        # Check lengths
        if not ((len(tasks) == len(functional)) or (len(functional) == 1)):
            raise ValueError('Functional must be assigned globally or per'
                             'task.')

        if not ((len(tasks) == len(basis_set)) or (len(basis_set) == 1)):
            raise ValueError('Basis set must be assigned globally or per'
                             'task.')

        if not ((len(tasks) == len(ao_basis)) or (len(ao_basis) == 1)):
            raise ValueError('AO basis must be assigned globally or per task.')

        if not ((len(tasks) == len(cosmo)) or (len(cosmo) == 1)):
            raise ValueError('Maximum iterations must be assigned globally or'
                             'per task.')

        # Extract atom index information
        self.idx = self.geom.get_atom_indices(atoms=atoms)

        # Generate header information
        config += self._configure_header(scratch_dir=scratch_dir,
                                         mem_global=mem_global,
                                         mem_heap=mem_heap,
                                         mem_stack=mem_stack)

        # Load geometry
        config += self._configure_load(charge=charge)

        # Configure tasks
        for task, f, b, a, c in zip(tasks, cycle(functional), cycle(basis_set),
                                    cycle(ao_basis), cycle(cosmo)):
            # TODO: finish this
            config += self.task_map[task](functional=f,
                                          basis_set=b,
                                          ao_basis=a,
                                          energy=energy,
                                          frequency=frequency,
                                          temp=temp,
                                          cosmo=c,
                                          solvent=solvent,
                                          gas=gas,
                                          max_iter=max_iter)

        # Store as atrribute
        self.config = config

        return self.config

    def configure_from_template(self, path, basename_override=None,
                                dirname_override=None, **kwargs):
        '''
        Configure NWChem simulation from template file.

        Use for NWChem functionality not exposed by the wrapper.
        Template contains ${keyword} entries that will be replaced by entries
        in `**kwargs`, if present. By default, ${basename} and ${dirname} must
        be included in the template and will be populated automatically.
        Override this behavior through use of appropriate keyword arguments.

        Parameters
        ----------
        path : str
            Path to template file.
        basename_override : str
            Override managed basename with user-supplied alternative.
        dirname_override : str
            Override managed directory name with user-supplied alternative.
        **kwargs
            Keyword arguments that will be subsituted in the template.

        Returns
        -------
        str
            NWChem configuration.

        '''

        # Add/override class-managed kwargs
        if basename_override is not None:
            kwargs['basename'] = basename_override
        else:
            kwargs['basename'] = self.geom.basename

        if dirname_override is not None:
            kwargs['dirname'] = dirname_override
        else:
            kwargs['dirname'] = self.temp_dir.name

        # Tied to save_geometry, required
        kwargs['fmt'] = self.fmt

        # Open template
        with open(path, 'r') as f:
            template = Template(f.read())

        # Store as attribute
        self.config = template.substitute(**kwargs)

        return self.config

    def save_config(self):
        '''
        Write generated NWChem configuration to file.

        '''

        # Write to file
        with open(os.path.join(self.temp_dir.name,
                               self.geom.basename + '.nw'), 'w') as f:
            f.write(self.config)

    def run(self):
        '''
        Run the NWChem simulation according to configured inputs.

        '''

        infile = os.path.join(self.temp_dir.name, self.geom.basename + '.nw')
        outfile = os.path.join(self.temp_dir.name, self.geom.basename + '.out')
        logfile = os.path.join(self.temp_dir.name, self.geom.basename + '.log')
        subprocess.call('nwchem {} > {} 2> {}'.format(infile,
                                                      outfile,
                                                      logfile), shell=True)

    def finish(self, keep_files=False, path=None):
        '''
        Parse NWChem simulation results and clean up temporary directory.

        Parameters
        ----------
        keep_files : bool
            Indicate whether to keep all intermediate files (relevant data will
            automatically be parsed).
        path : str
            Directory to copy intermediate files. Only used if `keep_files` is
            True.

        Returns
        -------
        :obj:`~isicle.parse.NWChemResult`
            Parsed result data.

        '''

        parser = NWChemParser()
        parser.load(os.path.join(self.temp_dir.name,
                                 self.geom.basename + '.out'))
        result = parser.parse(to_parse=['energy', 'shielding', 'spin',
                                        'molden', 'frequency'])

        if keep_files is True:
            import shutil
            import glob

            if path is None:
                raise ValueError('Must supply `path`.')
            else:
                # TODO: anything else to keep?
                shutil.copy2(os.path.join(self.temp_dir.name,
                                          self.geom.basename + '.nw'), path)
                shutil.copy2(os.path.join(self.temp_dir.name,
                                          self.geom.basename + '.out'), path)
                shutil.copy2(os.path.join(self.temp_dir.name,
                                          self.geom.basename + '.log'), path)

                geoms = glob.glob(os.path.join(self.temp_dir.name,
                                               '*.{}'.format(self.fmt)))

                [shutil.copy2(x, path) for x in geoms]

        # Remove temporary files
        self.temp_dir.cleanup()

        return result
