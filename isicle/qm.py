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
        :class:`~isicle.interfaces.WrapperInterface`.

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
    :obj:`~isicle.qm.QMWrapper`
        Instance of selected QMWrapper.

    '''

    # Select program
    return _program_selector(program).run(geom, template=template, **kwargs)


class NWChemWrapper(WrapperInterface):
    '''
    Wrapper for NWChem functionality.

    Implements :class:`~isicle.interfaces.WrapperInterface` to ensure
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

    def __init__(self, temp_dir=None):
        '''
        Initialize :obj:`~isicle.qm.NWChemWrapper` instance.

        Establishes aliases for preconfigured tasks.

        '''

        self.task_map = {'optimize': self._configure_optimize,
                         'energy': self._configure_energy,
                         'frequency': self._configure_frequency,
                         'shielding': self._configure_shielding,
                         'spin': self._configure_spin}

        self.task_order = {'optimize': 0,
                           'energy': 1,
                           'frequency': 2,
                           'shielding': 3,
                           'spin': 4}

        # Set up temporary directory
        self.temp_dir = tempfile.TemporaryDirectory(dir=temp_dir)


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

        # Save
        self.save_geometry()

    def save_geometry(self):
        '''
        Save internal :obj:`~isicle.geometry.Geometry` representation to file.

        Creates temporary directory for intermediate files.

        '''

        # Path operations
        outfile = os.path.join(self.temp_dir.name,
                               '{}.xyz'.format(self.geom.basename))

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
             'dirname': self.temp_dir.name,
             'charge': charge}

        return ('\ncharge {charge}\n'
                'geometry noautoz noautosym\n'
                ' load {dirname}/{basename}.xyz\n'
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
             'max_iter': max_iter}

        return ('\ndriver\n'
                ' maxiter {max_iter}\n'
                ' xyz {basename}_geom\n'
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

    def _configure_frequency(self, temp=298.15, basis_set='6-31G*', ao_basis='cartesian',
                             functional='b3lyp', cosmo=False, solvent='H2O', gas=False,
                             **kwargs):
        '''
        Configure frequency block of NWChem configuration.

        Includes basis and DFT blocks; can include COSMO block.

        Parameters
        ----------
        temp : float
            Temperature for frequency calculation.
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
        **kwargs
            Arbitrary additional arguments (unused).

        Returns
        -------
        str
            Frequency block of NWChem configuration.

        '''

         # Add basis block
        s = self._configure_basis(basis_set=basis_set, ao_basis=ao_basis)

        # Add DFT block
        s += self._configure_dft(functional=functional)

        # Add COSMO block
        if cosmo:
            s += self._configure_cosmo(solvent=solvent, gas=gas)

        # Add frequency block
        s += ('\nfreq\n'
              ' temp 1 {}\n'
              'end\n').format(temp)
        s += '\ntask dft freq ignore\n'

        return s

    def _configure_energy(self, basis_set='6-31G*', ao_basis='cartesian',
                          functional='b3lyp', cosmo=False, solvent='H2O', gas=False,
                          **kwargs):
        '''
        Configure energy block of NWChem configuration.

        Includes basis and DFT blocks; can include COSMO block.

        Parameters
        ----------
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
        **kwargs
            Arbitrary additional arguments (unused).

        Returns
        -------
        str
            Energy block of NWChem configuration.
        '''

         # Add basis block
        s = self._configure_basis(basis_set=basis_set, ao_basis=ao_basis)

        # Add DFT block
        s += self._configure_dft(functional=functional)

        # Add COSMO block
        if cosmo:
            s += self._configure_cosmo(solvent=solvent, gas=gas)

        s += '\ntask dft energy ignore\n'

        return s

    def _configure_optimize(self, basis_set='6-31G*', ao_basis='cartesian',
                            functional='b3lyp', max_iter=150,
                            cosmo=False, solvent='H2O', gas=False,
                            **kwargs):
        '''
        Generate meta optimization block of NWChem configuration.

        Includes basis, DFT, and driver blocks; can include COSMO block.

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
        s += self._configure_driver(max_iter=max_iter)

        # Add COSMO block
        if cosmo:
            s += self._configure_cosmo(solvent=solvent, gas=gas)

        # Add optimize task
        s += '\ntask dft optimize ignore\n'

        return s

    def _configure_shielding(self, basis_set='6-31G*', ao_basis='cartesian',
                             functional='b3lyp', cosmo=True, solvent='H2O',
                             gas=False, **kwargs):
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

        # Add property task
        s += '\ntask dft property ignore\n'

        return s

    def _configure_spin(self, max_pairs=30, basis_set='6-31G*',
                        ao_basis='cartesian', functional='b3lyp', cosmo=True,
                        solvent='H2O', gas=False, **kwargs):
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
                    s += '\ntask dft property ignore\n\n'
                s += '\nproperty\n'
                npairs = min(len(pairs) - i, max_pairs)
                s += ' SPINSPIN {}'.format(npairs)
            s += ' {} {}'.format(*p)

        s += '\nend\n'

        # Add property task
        s += '\ntask dft property ignore\n'

        return s

    def configure(self, tasks='energy', functional='b3lyp',
                  basis_set='6-31g*', ao_basis='cartesian', charge=0,
                  atoms=['C', 'H'], temp=298.15, cosmo=False, solvent='H2O',
                  gas=False, max_iter=150, mem_global=1600, mem_heap=100,
                  mem_stack=600, scratch_dir='/scratch'):
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
            Atom types of interest. Only used for `spin` and `shielding` tasks.
        temp : float
            Temperature for frequency calculation. Only used if `frequency` is
            a selected task.
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
        idx = self.geom.get_atom_indices(atoms=atoms)
        new_idx = []
        for i in idx:
            new_idx.append(str(int(i)+1))
        
        self.idx = new_idx

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
                                          temp=temp,
                                          cosmo=c,
                                          solvent=solvent,
                                          gas=gas,
                                          max_iter=max_iter)

        # Store as atrribute
        self.config = config

        # Save
        self.save_config()

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

        # Open template
        with open(path, 'r') as f:
            template = Template(f.read())

        # Store as attribute
        self.config = template.substitute(**kwargs)
        
        # Save
        self.save_config()

        return self.config

    def save_config(self):
        '''
        Write generated NWChem configuration to file.

        '''

        # Write to file
        with open(os.path.join(self.temp_dir.name,
                               self.geom.basename + '.nw'), 'w') as f:
            f.write(self.config)

    def submit(self):
        '''
        Submit the NWChem simulation according to configured inputs.

        '''

        infile = os.path.join(self.temp_dir.name, self.geom.basename + '.nw')
        outfile = os.path.join(self.temp_dir.name, self.geom.basename + '.out')
        logfile = os.path.join(self.temp_dir.name, self.geom.basename + '.log')
        subprocess.call('nwchem {} > {} 2> {}'.format(infile,
                                                      outfile,
                                                      logfile), shell=True)

    def finish(self):
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
        result = parser.parse()

        self.__dict__.update(result)
        self.geom.add_global_properties({k: v for k, v in result.items() if k != 'geom'})
        return self

    def run(self, geom, template=None, **kwargs):
        '''
        Optimize geometry via density functional theory using supplied functional
        and basis set.

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.
        template : str
            Path to optional template to bypass default configuration process.
        **kwargs
            Keyword arguments to configure the simulation.
            See :meth:`~isicle.qm.NWChemWrapper.configure`.

        Returns
        -------
        :obj:`~isicle.qm.NWChemWrapper`
            Wrapper object containing relevant outputs from the simulation.

        '''

        # Set geometry
        self.set_geometry(geom)

        # Configure
        if template is not None:
            self.configure_from_template(template)
        else:
            self.configure(**kwargs)

        # Run QM simulation
        self.submit()

        # Finish/clean up
        self.finish()

        return self
