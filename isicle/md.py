import os
import subprocess

import isicle
from isicle.geometry import Geometry, XYZGeometry
from isicle.interfaces import WrapperInterface
from isicle.parse import XTBParser

'''
Files resulting from an xtb job always run in the same directory that the command is
issued in, no matter where the input is. Can direct the .log file, but no other files.
'''


def _program_selector(program):
    '''
    Selects a supported molecular dynamics program for associated simulation.
    Currently only NWChem has been implemented.

    Parameters
    ----------
    program : str
        Alias for program selection (xtb).

    Returns
    -------
    program
        Wrapped functionality of the selected program. Must implement
        :class:`~isicle.interfaces.MDWrapperInterface`.

    '''

    program_map = {'xtb': XTBWrapper, 'rdkit': RDKitWrapper}

    if program.lower() in program_map.keys():
        return program_map[program.lower()]()
    else:
        raise ValueError(
            '{} not a supported molecular dynamics program.'.format(program))


def md(geom, program='xtb', **kwargs):
    '''
    Optimize geometry via molecular dyanmics using supplied forcefield
    and basis set.

    Parameters
    ----------
    program : str
        Alias for program selection (xtb).

    Returns
    -------
    result
        Object containing relevant outputs from the simulation.

    '''

    # Select program
    return _program_selector(program).run(geom, **kwargs)


class XTBWrapper(XYZGeometry, WrapperInterface):
    '''
    Wrapper for xtb functionality.

    Implements :class:`~isicle.interfaces.WrapperInterface` to ensure required methods are exposed.

    Attributes
    ----------
    temp_dir : str
        Path to temporary directory used for simulation.
    task_map : dict
        Alias mapper for supported molecular dynamic presets. Includes
        "optimize", "crest", "nmr", "protonate", "deprotonate", and "tautomer".
    geom : :obj:`isicle.geometry.Geometry`
        Internal molecule representation.
    fmt : str
        File extension indicator.
    job_list : str
        List of commands for simulation.

    '''

    _defaults = ['geom']
    _default_value = None

    def __init__(self, **kwargs):
        self.__dict__.update(dict.fromkeys(
            self._defaults, self._default_value))
        self.__dict__.update(**kwargs)

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
        self.basename = self.geom.basename

        # Save geometry
        self.save_geometry()

    def save_geometry(self, fmt='xyz'):
        '''
        Save internal :obj:`~isicle.geometry.Geometry` representation to file.

        Parameters
        ----------
        fmt : str
            Filetype used by xtb. Must be "xyz", "smi", ".inchi", ".mol", ".xyz",
            ".pdb", ".pkl".

        '''
        # Path operationspyth
        self.temp_dir = isicle.utils.mkdtemp()
        self.fmt = fmt.lower()
        geomfile = os.path.join(self.temp_dir,
                                '{}.{}'.format(self.basename,
                                               self.fmt.lower()))

        # All other formats
        isicle.io.save(geomfile, self.geom)
        self.geom.path = geomfile

    def _configure_xtb(self, forcefield='gfn2', optlevel='normal', charge=None, solvation=None):
        '''
        Set command line for xtb simulations.

        Parameters
        ----------
        forcefield : str
            GFN forcefield for the optimization
            Default: gff
            Supported forcefields: gfn2, gfn1, gff
        optlevel : str
            Optimization convergence level
            Default : normal
            Supported : crude, sloppy, loose, lax, normal, tight, vtight extreme
        charge : int
            Charge of molecular system.
            Default : 0 (Neutral charge)

        '''

        # Add base command
        s = 'xtb '

        # Add geometry
        s += '{}.{}'.format(self.basename, self.fmt.lower())

        # Add optimize tag
        s += ' --opt ' + optlevel + ' '

        # Add forcefield
        s += '--' + forcefield + ' '

        # Add optional charge
        if charge is not None:
            s += '--chrg ' + charge + ' '

        # Add optional implicit solvation
        if solvation is not None:
            s += '--alpb ' + solvation + ' '

        # Add output
        s += '&>' + ' '

        s += '{}.{}'.format(self.basename, "out")
        return s

    def _configure_crest(self, ewin=6, optlevel='Normal', forcefield='gfn2',
                         protonate=False, deprotonate=False, tautomerize=False,
                         ion=None, charge=None, dryrun=False, processes=1,
                         solvation=None, ignore_topology=False):
        '''
        Set command line for crest simulations.

        Parameters
        ----------
        ewin : int
            Energy window (kcal/mol) for conformer, (de)protomer, or tautomer search.
            Default : 6
        optlevel : str
            Optimization convergence level
            Default : normal
            Supported : crude, sloppy, loose, lax, normal, tight, vtight extreme
        forcefield : str
            GFN forcefield for the optimization
            Default: gff
            Supported forcefields: gfn2, gfn1, gff
        protonate : bool
            Signal to initiate protomer search. Suggested ewin = 30.
            Default : False
        deprotonate : bool
            Signal to initiate deprotonated conformers. Suggesting ewin = 30.
            Default : False
        tautomer : bool
            Signal to initiate tautomer search.
            Default : False
        ion : str
            Keyword to couple with protonate to ionize molecule with an ion other than a proton.
            See :obj:`~isicle.adduct.parse_ion` for list of ion options.
        charge : int
            Charge of molecular system.
            Default : 0 (Neutral charge)
        '''

        # Start base command
        s = 'crest '

        # Add geometry
        s += str(os.path.join(self.temp_dir,
                              '{}.{}'.format(self.basename,
                                             self.fmt.lower())))

        s += ' '
        # Add optional tag
        if protonate:
            s += '-protonate '
        elif deprotonate:
            s += '-deprotonate '
        elif tautomerize:
            s += '-tautomerize '

        if ion is not None:
            s += '-swel ' + ion + ' '

        if charge is not None:
            s += '-chrg ' + str(charge) + ' '

        # Add dryrun option
        if dryrun:
            s += '--dryrun '

        # Add energy window
        s += '--ewin ' + str(ewin) + ' '

        # Add optlevel
        s += '--optlevel ' + optlevel + ' '

        # Add forcefield
        s += '-' + forcefield + ' '

        # Add optional solvation
        if solvation is not None:
            s += '--alpb ' + solvation + ' '

        # Number of processes
        s += '-T ' + str(processes) + ' '

        if ignore_topology:
            s += '--noreftopo '

        # Add output
        s += '&>' + ' '

        s += os.path.join(self.temp_dir,
                          '{}.{}'.format(self.basename,
                                         "out"))

        return s

    def configure(self, task='optimize', forcefield='gfn2', charge=None,
                  ewin=6, ion=None, optlevel='Normal', dryrun=False, processes=1,
                  solvation=None, ignore_topology=False):
        '''
        Generate command line

        Parameters
        ----------
        tasks : str
            Set task to "optimize", "conformer", "protonate", "deprotonate", or "tautomerize".
            Default : "optimize"
        forcefield : str
            GFN forcefield for the optimization
            Default: gff
            Supported forcefields: gfn2, gfn1, gff
        ewin : int
            Energy window (kcal/mol) for conformer(set to 6), (de)protomer(set to 30), or tautomer(set to 30) search.
            Default : 6
        ion : str
            Ion for protomer calculation.
        optlevel : str or list of str
            Set optimization level. Supply globally or per task.
        ion : str
            Keyword to couple with protonate to ionize molecule with an ion other than a proton.
            See :obj:`~isicle.adduct.parse_ion` for list of ion options.
        charge : int
            Charge of molecular system.
            Default : 0 (Neutral charge)
        '''

        if type(task) == list:
            raise TypeError('Initiate one xtb or crest job at a time.')
        if type(forcefield) == list:
            raise TypeError('Initiate one forcefield at a time.')
        if type(optlevel) == list:
            raise TypeError('Initiate one opt level at a time.')

        if task == 'optimize':
            config = self._configure_xtb(optlevel=optlevel,
                                         forcefield=forcefield)

        else:
            if task == 'conformer':
                p, d, t, i = False, False, False, None

            elif task == 'protonate':
                p, d, t, i = True, False, False, ion

            elif task == 'deprotonate':
                p, d, t, i = False, True, False, ion

            elif task == 'tautomerize':
                p, d, t, i = False, False, True, ion

            if p is not None:
                config = self._configure_crest(ewin=ewin,
                                               optlevel=optlevel,
                                               forcefield=forcefield,
                                               protonate=p,
                                               deprotonate=d,
                                               tautomerize=t,
                                               ion=i,
                                               charge=charge,
                                               dryrun=dryrun,
                                               processes=processes,
                                               solvation=solvation,
                                               ignore_topology=ignore_topology)
            else:
                raise Error(
                    'Task not assigned properly, please choose optimize, conformer, protonate, deprotonate, or tautomerize')

        self.task = task

        self.config = config

    def submit(self):
        '''
        Run xtb or crest simulation according to configured inputs.
        '''
        owd = os.getcwd()
        os.chdir(self.temp_dir)
        job = self.config
        subprocess.call(job, shell=True)
        os.chdir(owd)

    def finish(self):
        '''
        Parse results, save xtb output files, and clean temporary directory
        '''

        parser = XTBParser()

        parser.load(os.path.join(self.temp_dir, self.basename + '.out'))
        self.output = parser.load(os.path.join(
            self.temp_dir, self.basename + '.out'))

        result = parser.parse()

        self.__dict__.update(result)

        for i in self.geom:
            i.add___dict__({k: v for k, v in result.items() if k != 'geom'})
            i.__dict__.update(basename=self.basename)

        if self.task != 'optimize':
            conformerID = 1
            for i in self.geom:
                i.__dict__.update(conformerID=conformerID)
                conformerID += 1
            return self
        else:
            self.geom = self.geom[0]

    def run(self, geom, **kwargs):
        '''
        Optimize geometry via density functional theory using supplied functional
        and basis set.

        Parameters
        ----------
        geom : :obj:`~isicle.geometry.Geometry`
            Molecule representation.
        **kwargs
            Keyword arguments to configure the simulation.
            See :meth:`~isicle.md.XTBWrapper.configure`.

        Returns
        -------
        :obj:`~isicle.md.XTBWrapper`
            Wrapper object containing relevant outputs from the simulation.

        '''

        # New instance
        self = XTBWrapper()

        # Set geometry
        self.set_geometry(geom)

        # Configure
        self.configure(**kwargs)

        # Run QM simulation
        self.submit()

        # Finish/clean up
        self.finish()

        return self

    def get_structures(self):
        '''
        Extract all structures from containing object as a conformational ensemble.

        Returns
        -------
        :obj:`~isicle.conformers.ConformationalEnsemble`
            Conformational ensemble.

        '''
        if isinstance(self.geom, isicle.conformers.ConformationalEnsemble):
            return self.geom

        raise TypeError(
            'Object does not contain multiple structures. Use `get_structure` instead.')

    def get_structure(self):
        '''
        Extract structure from containing object.

        Returns
        -------
        :obj:`~isicle.geometry.XYZGeometry`
            Structure instance. 

        '''
        if isinstance(self.geom, isicle.conformers.ConformationalEnsemble):
            raise TypeError(
                'Object contains multiple structures. Use `get_structures` instead.')

        return self.geom


class RDKitWrapper(Geometry, WrapperInterface):

    '''
    Wrapper for rdkit functionality.

    Implements :class:`~isicle.interfaces.WrapperInterface` to ensure required methods are exposed.

    Attributes
    ----------
    temp_dir : str
        Path to temporary directory used for simulation.
    task_map : dict
        Alias mapper for supported molecular dynamic presets. Includes
        "optimize", "conformer", "nmr", "protonate", "deprotonate", and "tautomer".
    geom : :obj:`isicle.geometry.Geometry`
        Internal molecule representation.
    fmt : str
        File extension indicator.
    job_list : str
        List of commands for simulation.

    '''

    _defaults = ['geom']
    _default_value = None

    def __init__(self, **kwargs):
        self.__dict__.update(dict.fromkeys(
            self._defaults, self._default_value))
        self.__dict__.update(**kwargs)

    def set_geometry(self):
        return

    def configure(self):
        return

    def submit(self):
        return

    def run(self):
        return

    def finish(self):
        return


def amber():
    # Don't work on this one yet
    # TODO: define input, arguments, and output
    # Should return instance (or list) of MDOptimizedGeometry
    raise NotImplementedError

# TODO: add other MD methods as needed
