from isicle.interfaces import WrapperInterface
from isicle.parse import XTBParser
import tempfile
import os

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

    program_map = {'xtb': XTBWrapper}

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
    geom : :obj:`~isicle.geometry.Geometry`
        Molecule representation.
    program : str
        Alias for program selection (xtb).
    **kwargs
        Keyword arguments to configure the simulation.
        See :meth:`~isicle.qm.XTBWrapper.configure`.

    Returns
    -------
    result
        Object containing relevant outputs from the simulation.

    '''

    # Select program
    mdw = _program_selector(program)

    # Set geometry
    mdw.set_geometry(geom)

    # Save geometry
    if 'fmt' in kwargs:
        mdw.save_geometry(fmt=kwargs.pop('fmt'))
    else:
        mdw.save_geometry()

    # Build command line 
    mdw.configure()

    #Save configuration
    mdw.save_config()

    # Run MD simulation
    mdw.run()

    # Create new Geometry with updated structure
    # res['geometry'] will be None or a path to an xyz file(s).
    res = mdw.finish()

    geom = geom._update_structure(False, xyz=res.geometry)

    # Erase old properties and add new event and MD properties
    geom.global_properties = {}
    geom._update_history('md')
    geom.add_global_properties(res.to_dict())

    # Finish/clean up
    return geom, res

class XTBWrapper(WrapperInterface):
    '''
    Wrapper for xtb functionality.

    Implements :class:`~isicle.interfaces.MDWrapperInterface` to ensure required methods are exposed.

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

    def __init__(self):
        '''
        Initialize :obj:`MDChemWrapper` instance.

        Creates temporary directory for intermediate files, establishes aliases
        for preconfigured tasks.

        '''
        self.temp_dir = tempfile.TemporaryDirectory()
        self.task_map = {'xtb': self._configure_xtb,
                         'crest': self._configure_crest}

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
            Filetype used by xtb. Must be "xyz", "smi", ".inchi", ".mol", ".xyz",
            ".pdb", ".pkl".

        '''
        # Path operationspyth
        self.fmt = fmt.lower()
        outfile = os.path.join(self.temp_dir.name,
                               '{}.{}'.format(self.geom.basename,
                                              self.fmt.lower()))

        # All other formats
        self.geom.save(outfile)
        self.geom.path = outfile

    def _configure_xtb(self, forcefield='gff', optlevel='normal',charge=None):
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
        s += '{}.{}'.format(self.geom.basename,self.fmt.lower())
 
        # Add optimize tag
        s += ' --opt ' + optlevel + ' '

        # Add forcefield
        s += '--'+ forcefield + ' ' 

        # Add optional charge
        if charge is not None:
            s += '--chrg '+ charge + ' '

        # Add output
        s += '&>' + ' ' 

        s += '{}.{}'.format(self.geom.basename, "out")
        return s


    def _configure_crest(self, ewin=6, optlevel='Normal', forcefield='gff',
                       protonate=False, deprotonate=False, tautomerize=False,
                       ion=None, charge=None, dryrun=False):

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
        s += str(os.path.join(self.temp_dir.name,
                               '{}.{}'.format(self.geom.basename,
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
            s += '-chrg ' + charge + ' '

        # Add dryrun option
        if dryrun:
            s += '--dryrun '

        # Add energy window
        s += '--ewin ' + str(ewin) + ' '

        # Add optlevel
        s += '--optlevel ' + optlevel + ' '

        # Add forcefield
        s += '-'+ forcefield + ' '

        # Add scratch folder
        s += '--scratch '

        # Add output
        s += '&>' + ' '

        s += os.path.join(self.temp_dir.name,
                               '{}.{}'.format(self.geom.basename,
                                              "out"))

        return s

    def configure(self, task='optimize', forcefield='gff', charge=None,
                  ewin=6, ion=None, optlevel='Normal', dryrun=False):
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

            config = self._configure_crest(ewin=ewin, 
                                           optlevel=optlevel, 
                                           forcefield=forcefield,
                                           protonate=p, 
                                           deprotonate=d, 
                                           tautomerize=t,
                                           ion=i,
                                           charge=charge,
                                           dryrun=dryrun)

        self.task = task

        self.config = config

    def save_config(self):
        '''Filler function to match WrapperInterface'''
        self = self
        return 

    def run(self):
        '''
        Run xtb or crest simulation according to configured inputs.
        '''
        owd = os.getcwd()
        os.chdir(self.temp_dir.name)
        job = self.config
        os.system(job)

        os.chdir(owd)

    def finish(self, keep_files=False, path=None):
        '''
        Parse results, save xtb output files, and clean temporary directory

        '''

        parser = XTBParser()

        parser.load(os.path.join(self.temp_dir.name, self.geom.basename + '.out'))
        result = parser.parse(to_parse=['energy', 'geometry'])

        if keep_files is True:
            import shutil

            if path is None:
                raise ValueError('Must supply `path`.')
            else:
                shutil.copy2(os.path.join(self.temp_dir.name,
                                          self.geom.basename + '.out'), path)
                if 'optimize' in self.task:
                    shutil.copy2(os.path.join(self.temp_dir.name,
                                              'xtbopt.{}'.format(self.fmt.lower())), path)
                if 'crest' in self.task:
                    shutil.copy2(os.path.join(self.temp_dir.name,
                                              'crest_conformers.xyz'), path)
                if 'protonate' in self.task:
                    shutil.copy2(os.path.join(self.temp_dir.name,
                                              'protonated.xyz'), path)
                if 'deprotonate' in self.task:
                    shutil.copy2(os.path.join(self.temp_dir.name,
                                              'deprotonated.xyz'), path)
                if 'tautomer' in self.task:
                    shutil.copy2(os.path.join(self.temp_dir.name,
                                              'tautomers.xyz'), path)

        # Remove temporary files
        self.temp_dir.cleanup()

        return result

def amber():
    # Don't work on this one yet
    # TODO: define input, arguments, and output
    # Should return instance (or list) of MDOptimizedGeometry
    raise NotImplementedError

# TODO: add other MD methods as needed
