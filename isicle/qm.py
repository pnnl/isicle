import glob
import os
import subprocess
from collections import OrderedDict
from itertools import cycle
from string import Template

import isicle
from isicle.interfaces import WrapperInterface
from isicle.utils import safelist


def _backend_selector(backend):
    """
    Selects a supported quantum mechanical backend for associated simulation.
    Currently only NWChem and ORCA have been implemented.

    Parameters
    ----------
    backend : str
        Alias for backend selection (e.g. NWChem, ORCA).

    Returns
    -------
    backend
        Wrapped functionality of the selected backend. Must implement
        :class:`~isicle.interfaces.WrapperInterface`.

    """

    backend_map = {'nwchem': NWChemWrapper,
                   'orca': ORCAWrapper}

    if backend.lower() in backend_map.keys():
        return backend_map[backend.lower()]()
    else:
        raise ValueError(('{} not a supported quantum mechanical backend.')
                         .format(backend))


def dft(geom, backend='NWChem', **kwargs):
    """
    Perform density functional theory calculations according to supplied task list
    and configuration parameters.

    Parameters
    ----------
    geom : :obj:`~isicle.geometry.Geometry`
        Molecule representation.
    backend : str
        Alias for backend selection (NWChem, ORCA).
    kwargs
        Keyword arguments passed to selected backend.

    Returns
    -------
    :obj:`~isicle.qm.QMWrapper`
        Wrapper object containing relevant outputs from the simulation.

    """

    # Select backend
    return _backend_selector(backend).run(geom, **kwargs)


class NWChemWrapper(WrapperInterface):
    """
    Wrapper for NWChem functionality.

    Attributes
    ----------
    temp_dir : str
        Path to temporary directory used for simulation.
    geom : :obj:`~isicle.geometry.Geometry`
        Internal molecule representation.
    result : dict
        Dictionary containing simulation results.

    """

    _defaults = ["geom", "result", "temp_dir"]
    _default_value = None

    def __init__(self):
        """
        Initialize :obj:`~isicle.qm.NWChemWrapper` instance.

        """

        self._task_map = {'optimize': self._configure_optimize,
                         'energy': self._configure_energy,
                         'frequency': self._configure_frequency,
                         'shielding': self._configure_shielding,
                         'spin': self._configure_spin}

        self._task_order = {'optimize': 0,
                           'energy': 1,
                           'frequency': 2,
                           'shielding': 3,
                           'spin': 4}

        # Set default attributes
        self.__dict__.update(dict.fromkeys(self._defaults, self._default_value))

        # Set up temporary directory
        self.temp_dir = isicle.utils.mkdtemp()

    def set_geometry(self, geom):
        """
        Set :obj:`~isicle.geometry.Geometry` instance for simulation.

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.

        """

        # Assign geometry
        self.geom = geom.__copy__()

        # Save
        self._save_geometry()

    def _save_geometry(self):
        """
        Save internal :obj:`~isicle.geometry.Geometry` representation to file.

        """

        # Path operations
        geomfile = os.path.join(self.temp_dir,
                                '{}.xyz'.format(self.geom.basename))

        # Save
        isicle.save(geomfile, self.geom)
        self.geom.path = geomfile

    def _configure_header(self, scratch_dir=None, mem_global=1600,
                          mem_heap=100, mem_stack=600):
        """
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

        """

        d = {'basename': self.geom.basename,
             'dirname': self.temp_dir,
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

    def _configure_load(self):
        """
        Generate geometry load block of NWChem configuration.

        Returns
        -------
        str
            Load block of NWChem configuration.

        """

        d = {'basename': self.geom.basename,
             'dirname': self.temp_dir,
             'charge': self.geom.formal_charge}

        return ('\ncharge {charge}\n'
                'geometry noautoz noautosym\n'
                ' load {dirname}/{basename}.xyz\n'
                'end\n').format(**d)

    def _configure_basis(self, basis_set='6-31G*', ao_basis='cartesian'):
        """
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

        """

        d = {'ao_basis': ao_basis,
             'basis_set': basis_set}

        return ('\nbasis {ao_basis}\n'
                ' * library {basis_set}\n'
                'end\n').format(**d)

    def _configure_dft(self, functional='b3lyp', odft=False):
        """
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

        """

        d = {'functional': functional,
             'dft': 'odft' if odft is True else 'dft'}

        s = '\ndft\n'

        if odft is True:
            s += ' odft\n'

        s += ' xc {functional}\n'.format(**d)
        s += ' mulliken\n'               # Do we need this line?
        s += ' print "mulliken ao"\n'    # (and this one?)
        s += 'end\n'

        return s

    def _configure_driver(self, max_iter=150):
        """
        Generate driver block of NWChem configuration.

        Parameters
        ----------
        max_iter : int
            Maximum number of optimization iterations.

        Returns
        -------
        str
            Driver block of NWChem configuration.

        """

        d = {'basename': self.geom.basename,
             'max_iter': max_iter}

        return ('\ndriver\n'
                ' maxiter {max_iter}\n'
                ' xyz {basename}_geom\n'
                'end\n').format(**d)

    def _configure_cosmo(self, solvent='H2O', gas=False):
        """
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

        """

        d = {'solvent': solvent,
             'gas': gas}

        return ('\ncosmo\n'
                ' do_gasphase {gas}\n'
                ' solvent {solvent}\n'
                'end\n').format(**d)

    def _configure_frequency(self, temp=298.15, basis_set='6-31G*', ao_basis='cartesian',
                             functional='b3lyp', cosmo=False, solvent='H2O', gas=False,
                             **kwargs):
        """
        Configure frequency block of NWChem configuration.

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

        Returns
        -------
        str
            Frequency block of NWChem configuration.

        """

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
        """
        Configure energy block of NWChem configuration.

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

        Returns
        -------
        str
            Energy block of NWChem configuration.

        """

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
                            cosmo=False, solvent='H2O', gas=False, **kwargs):
        """
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
            Indicate whether to u se gas phase calculations. Only used if
            `cosmo` is True.

        Returns
        -------
        str
            Optimization meta block of NWChem configuration.

        """

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
        """
        Generate meta shielding block of NWChem configuration.

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

        Returns
        -------
        str
            Shielding meta block of NWChem configuration.

        """

        # Extract atom index information
        # idx = self.geom.mol.get_atom_indices(atoms=atoms)
        # new_idx = []
        # for i in idx:
        #    new_idx.append(str(int(i)+1))

        # self.idx = new_idx

        # Add basis block
        s = self._configure_basis(basis_set=basis_set, ao_basis=ao_basis)

        # Add DFT block
        s += self._configure_dft(functional=functional)

        # Add COSMO block
        if cosmo:
            s += self._configure_cosmo(solvent=solvent, gas=gas)

        s += ('\nproperty\n'
              ' SHIELDING\n'
              'end\n')

        # Add property task
        s += '\ntask dft property ignore\n'

        return s

    def _configure_spin(self, bonds=1, basis_set='6-31G*',
                        ao_basis='spherical', functional='b3lyp', cosmo=True,
                        solvent='H2O', gas=False, **kwargs):
        """
        Generate meta spin-spin coupling block of NWChem configuration.

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

        Returns
        -------
        str
            Spin-spin coupling meta block of NWChem configuration.

        """

        def generate_pairs(mol, bonds=bonds):

            from rdkit import Chem

            matrix = Chem.GetAdjacencyMatrix(mol)

            pair_list = []

            for idx, atom in enumerate(matrix):
                # First round of neighbors
                neighbors = [i for i, x in enumerate(atom) if x == 1]

                for n in neighbors:
                    # Check if pair in pair_list, add one to index to match NWChem numbering.
                    if [idx+1, n+1] not in pair_list and [n+1, idx+1] not in pair_list:
                        pair_list.append([idx+1, n+1])

                    # Second round of neighbors
                    n_neighbors = [
                        i for i, x in enumerate(matrix[n]) if x == 1]

                    if bonds >= 2:
                        for n_n in n_neighbors:

                            if n_n != idx and atom[n_n] == 0:
                                if [idx+1, n_n+1] not in pair_list and [n_n+1, idx+1] not in pair_list:
                                    pair_list.append([idx+1, n_n+1])
                                atom[n_n] += 2

                            # Third line of neighbors
                            nn_neighbors = [i for i, x in enumerate(
                                matrix[n_n]) if x == 1]

                            if bonds >= 3:
                                for nn_n in nn_neighbors:
                                    if nn_n != idx and atom[nn_n] == 0:
                                        if [idx+1, nn_n+1] not in pair_list and [nn_n+1, idx+1] not in pair_list:
                                            pair_list.append([idx+1, nn_n+1])
                                        atom[nn_n] += 3
                            else:
                                continue
                        else:
                            continue

            # Build string of pairs for input
            s = ""
            for pair in pair_list:
                s += str(pair[0]) + " " + str(pair[1]) + " "

            return len(pair_list), s

        pair_count, pairs = generate_pairs(self.geom.mol, bonds=bonds)

        d = {'pair_count': pair_count, 'pairs': pairs}

        # Add basis block
        s = self._configure_basis(basis_set=basis_set, ao_basis=ao_basis)

        # Add DFT block
        s += self._configure_dft(functional=functional, odft=True)

        # Add COSMO block
        if cosmo:
            s += self._configure_cosmo(solvent=solvent, gas=gas)

        # Add spin block
        s += '\nend\n'
        s += '\ntask dft ignore\n\n'
        s += '\nproperty\n'
        s += ' SPINSPIN {pair_count} {pairs}'.format(**d)
        s += '\nend\n'

        # Add property task
        s += '\ntask dft property ignore\n'

        return s

    def configure(self, tasks='energy', functional='b3lyp',
                  basis_set='6-31g*', ao_basis='cartesian',
                  atoms=['C', 'H'], bonds=1, temp=298.15, cosmo=False, solvent='H2O',
                  gas=False, max_iter=150, mem_global=1600, mem_heap=100,
                  mem_stack=600, scratch_dir=None, processes=12):
        """
        Configure NWChem simulation.

        Parameters
        ----------
        tasks : str or list of str
            Tasks text.
            Tasks text.
        functional : str or list of str
            Functional selection. Supply globally or per task.
        basis_set : str or list of str
            Basis set selection. Supply globally or per task.
        ao_basis : str or list of str
            Angular function selection ("spherical", "cartesian"). Supply
            globally or per task.
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

        """

        # Cast to list safely
        tasks = safelist(tasks)
        functional = safelist(functional)
        basis_set = safelist(basis_set)
        ao_basis = safelist(ao_basis)
        atoms = safelist(atoms)
        cosmo = safelist(cosmo)
        solvent = safelist(solvent)

        # Set scratch directory
        if scratch_dir is None:
            scratch_dir = isicle.utils.mkdtemp()

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

        if not ((len(tasks) == len(solvent)) or (len(solvent) == 1)):
            raise ValueError('Solvents must be assigned globally or'
                             'per task.')

        # Generate header information
        config += self._configure_header(scratch_dir=scratch_dir,
                                         mem_global=mem_global,
                                         mem_heap=mem_heap,
                                         mem_stack=mem_stack)

        # Load geometry
        config += self._configure_load()

        # Configure tasks
        for task, f, b, a, c, so in zip(tasks, cycle(functional), cycle(basis_set),
                                        cycle(ao_basis), cycle(cosmo), cycle(solvent)):
            # TODO: finish this
            config += self._task_map[task](
                functional=f,
                basis_set=b,
                ao_basis=a,
                temp=temp,
                cosmo=c,
                gas=gas,
                max_iter=max_iter,
                solvent=so,
                bonds=bonds
                )

        # Store tasks as attribute
        self._tasks = tasks

        # Store number of processes as attribute
        self._processes = processes

        # Store as atrribute
        self._config = config

        # Save
        self.save_config()

        return self._config

    def configure_from_template(self, path, basename_override=None,
                                dirname_override=None, **kwargs):
        """
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

        """

        # Add/override class-managed kwargs
        if basename_override is not None:
            kwargs['basename'] = basename_override
        else:
            kwargs['basename'] = self.geom.basename

        if dirname_override is not None:
            kwargs['dirname'] = dirname_override
        else:
            kwargs['dirname'] = self.temp_dir

        # Open template
        with open(path, 'r') as f:
            template = Template(f.read())

        # Store as attribute
        self.config = template.substitute(**kwargs)

        # Save
        self.save_config()

        return self.config

    def save_config(self):
        """
        Write generated NWChem configuration to file.

        """

        # Write to file
        with open(os.path.join(self.temp_dir,
                               self.geom.basename + '.nw'), 'w') as f:
            f.write(self._config)

    def submit(self):
        """
        Submit the NWChem simulation according to configured inputs.

        """

        infile = os.path.join(self.temp_dir, self.geom.basename + '.nw')
        outfile = os.path.join(self.temp_dir, self.geom.basename + '.out')
        logfile = os.path.join(self.temp_dir, self.geom.basename + '.log')

        s = 'mpirun -n {} nwchem {} > {} 2> {}'.format(self._processes,
                                                       infile,
                                                       outfile,
                                                       logfile)

        subprocess.call(s, shell=True)

    def finish(self):
        """
        Parse NWChem simulation results.

        Returns
        -------
        dict
            Dictionary containing simulation results.

        """

        # Get list of outputs
        outfiles = glob.glob(os.path.join(self.temp_dir, '*'))

        # Result container
        result = {}

        # Split out geometry files
        geomfiles = sorted([x for x in outfiles if x.endswith('.xyz')])
        outfiles = sorted([x for x in outfiles if not x.endswith('.xyz')])

        # Enumerate geometry files
        result['xyz'] = OrderedDict()
        for geomfile in geomfiles:
            geom = isicle.load(geomfile)
            
            if '_geom-' in geomfile:
                idx = int(os.path.basename(geomfile).split('-')[-1].split('.')[0])
                result['xyz'][idx] = geom

            else:
                result['xyz']['input'] = geom

        # Rename final geometry
        result['xyz']['final'] = list(result['xyz'].values())[-1]

        # Enumerate output files
        for outfile in outfiles:
            # Split name and extension
            basename, ext = os.path.basename(outfile).rsplit('.', 1)

            # Read output content
            with open(outfile, 'rb') as f:
                contents = f.read()

            # Attempt utf-8 decode
            try:
                result[ext] = contents.decode('utf-8')
            except UnicodeDecodeError:
                result[ext] = contents

        # Assign to attribute
        self.result = result
        return self.result

    def parse(self):
        """
        Parse NWChem simulation results.

        Returns
        -------
        dict
            Dictionary containing parsed outputs from the simulation.

        """

        if self.result is None:
            raise RuntimeError("Must complete NWChem simulation.")

        parser = isicle.parse.NWChemParser(data=self.result)

        return parser.parse()

    def run(self, geom, template=None, tasks='energy', functional='b3lyp',
            basis_set='6-31g*', ao_basis='cartesian',
            atoms=['C', 'H'], bonds=1, temp=298.15, cosmo=False, solvent='H2O',
            gas=False, max_iter=150, mem_global=1600, mem_heap=100,
            mem_stack=600, scratch_dir=None, processes=12):
        """
        Perform density functional theory calculations according to supplied task list
        and configuration parameters.

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.
        template : str
            Path to optional template to bypass default configuration process.
        tasks : str or list of str
            List of calculations to perform. One or more of "optimize", "energy",
            "frequency", "shielding", "spin".
        functional : str or list of str
            Functional selection. Supply globally or per task.
        basis_set : str or list of str
            Basis set selection. Supply globally or per task.
        ao_basis : str or list of str
            Angular function selection ("spherical", "cartesian"). Supply
            globally or per task.
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
        dict
            Dictionary containing simulation results.

        """

        # Set geometry
        self.set_geometry(geom)

        # Configure
        if template is not None:
            self.configure_from_template(template)
        else:
            self.configure(tasks=tasks,
                           functional=functional, basis_set=basis_set,
                           ao_basis=ao_basis,atoms=atoms, bonds=bonds,
                           temp=temp, cosmo=cosmo, solvent=solvent, gas=gas,
                           max_iter=max_iter, mem_global=mem_global,
                           mem_heap=mem_heap, mem_stack=mem_stack,
                           scratch_dir=scratch_dir, processes=processes)

        # Run QM simulation
        self.submit()

        # Finish/clean up
        self.finish()

        return self

    def save(self, path):
        """
        Save result as pickle file.

        Parameters
        ----------
        path : str
            Path to output file.

        """

        if hasattr(self, 'result'):
            isicle.io.save_pickle(path, self.result)
        else:
            raise AttributeError("Object must have `result` attribute")


class ORCAWrapper(WrapperInterface):
    """
    Wrapper for ORCA functionality.

    Implements :class:`~isicle.interfaces.WrapperInterface` to ensure
    required methods are exposed.

    Attributes
    ----------
    temp_dir : str
        Path to temporary directory used for simulation.
    geom : :obj:`~isicle.geometry.Geometry`
        Internal molecule representation.
    config : str
        Configuration information for simulation.
    result : dict
        Dictionary containing simulation results.

    """

    _defaults = ["geom", "result", "temp_dir"]
    _default_value = None

    def __init__(self):
        """
        Initialize :obj:`~isicle.qm.ORCAWrapper` instance.

        """
        
        # Set defaults
        self.__dict__.update(dict.fromkeys(self._defaults, self._default_value))

        # Set up temporary directory
        self.temp_dir = isicle.utils.mkdtemp()

    def set_geometry(self, geom):
        """
        Set :obj:`~isicle.geometry.Geometry` instance for simulation.

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.

        """

        # Assign geometry
        self.geom = geom.__copy__()

        # Save
        self._save_geometry()

    def _save_geometry(self):
        """
        Save internal :obj:`~isicle.geometry.Geometry` representation to file in
        XYZ format.

        """

        # Path to output
        geomfile = os.path.join(self.temp_dir,
                                '{}.xyz'.format(self.geom.basename))

        # Save
        isicle.save(geomfile, self.geom)

        # Store path
        self.geom.path = geomfile

    def configure(self, simple_input=[], block_input={}, spin_multiplicity=1, processes=1, **kwargs):
        """
        Configure ORCA simulation.

        Parameters
        ----------
        simple_input : list
            List of simple input keywords. See `this <https://sites.google.com/site/orcainputlibrary/general-input>`__
            section of the ORCA docs.
        block_input : dict
            Dictionary defining configuration "blocks". Use names of blocks as keys, 
            lists of each block's content as values. To configure a line of block content
            directly, include as a complete string. Include key:value pairs as tuples. 
            See `this <https://sites.google.com/site/orcainputlibrary/general-input>`__
            section of the ORCA docs.
        spin_multiplicity : int
            Spin multiplicity of the molecule.
        processes : int
            Number of parallel processes.
        kwargs
            Additional keyword arguments fed as key:value pairs.

        Returns
        -------
        str
            ORCA configuration text.

        """

        # Safely cast to list
        simple_input = isicle.utils.safelist(simple_input)

        # Expand simple inputs
        config = '! ' + ' '.join(simple_input) + '\n'

        # Add processes
        if processes > 1:
            config += '%PAL NPROCS {} END\n'.format(processes)

        # Add geometry context
        config += '* xyzfile {:d} {:d} {}\n'.format(self.geom.formal_charge, spin_multiplicity, self.geom.path)

        # Expand keyword args
        for k, v in kwargs.items():
            config += '%{} {}\n'.format(k, v)

        # Expand block inputs
        for block, params in block_input.items():
            # Safely cast to list
            params = isicle.utils.safelist(params)

            # Block header
            block_text = '%{}\n'.format(block)

            # Block configuration
            for param in params:
                if type(param) is str:
                    block_text += param + '\n'
                elif type(param) is tuple:
                    block_text += ' '.join(map(str, param)) + '\n'
                else:
                    raise TypeError

            # End block
            block_text += 'end\n'

            # Append block to config
            config += block_text

        # Assign to attribute
        self.config = config

        # Save configuration
        self._save_config()

        return self.config

    def _save_config(self):
        """
        Write generated ORCA configuration to file.

        """

        # Write to file
        with open(os.path.join(self.temp_dir,
                               self.geom.basename + '.inp'), 'w') as f:
            f.write(self.config)

    def submit(self):
        """
        Submit the ORCA simulation according to configured inputs.

        """

        infile = os.path.join(self.temp_dir, self.geom.basename + '.inp')
        outfile = os.path.join(self.temp_dir, self.geom.basename + '.out')
        logfile = os.path.join(self.temp_dir, self.geom.basename + '.log')

        s = '`which orca` {} > {} 2> {}'.format(infile,
                                                outfile,
                                                logfile)

        subprocess.call(s, shell=True)

    def finish(self):
        """
        Collect ORCA simulation results.

        Returns
        -------
        dict
            Dictionary containing relevant outputs from the simulation.

        """

        # Get list of outputs
        outfiles = glob.glob(os.path.join(self.temp_dir, '*'))

        # Filter out temp files
        outfiles = [x for x in outfiles if not x.endswith('.tmp')]

        # Result container
        result = {}

        # Enumerate output files
        for outfile in outfiles:
            # Split name and extension
            basename, ext = os.path.basename(outfile).rsplit('.', 1)

            # Grab suitable variable names
            if any(basename.endswith(x) for x in ['_property', '_trj']):
                var_name = basename.split('_')[-1]
            else:
                var_name = ext

            # Load geometry
            if var_name == 'xyz':
                # Load xyz geometry
                geom = isicle.load(outfile)
                
                # Update coordinates of starting geometry
                geom = self.geom.update_coordinates(geom)

                # Add to result object
                result[var_name] = geom

            # Load other files
            else:
                # Read output content
                with open(outfile, 'rb') as f:
                    contents = f.read()

                # Attempt utf-8 decode
                try:
                    result[var_name] = contents.decode('utf-8')
                except UnicodeDecodeError:
                    result[var_name] = contents

        # Assign to attribute
        self.result = result
        return self.result

    def parse(self):
        """
        Parse ORCA simulation results.

        Returns
        -------
        dict
            Dictionary containing parsed outputs from the simulation.

        """

        if self.result is None:
            raise RuntimeError("Must complete ORCA simulation.")

        parser = isicle.parse.ORCAParser(data=self.result)

        return parser.parse()

    def run(self, geom, **kwargs):
        """
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
            See :meth:`~isicle.qm.ORCAWrapper.configure`.

        Returns
        -------
        dict
            Dictionary containing relevant outputs from the simulation.

        """

        # Set geometry
        self.set_geometry(geom)

        # Configure
        self.configure(**kwargs)

        # Run QM simulation
        self.submit()

        # Finish/clean-up outputs
        self.finish()

        return self

    def save(self, path):
        """
        Save result as pickle file.

        Parameters
        ----------
        path : str
            Path to output file.

        """

        if hasattr(self, 'result'):
            isicle.io.save_pickle(path, self.result)
        else:
            raise AttributeError("Object must have `result` attribute")
