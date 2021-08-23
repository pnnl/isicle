from numpy.core.records import array
from isicle.interfaces import WrapperInterface
from isicle.geometry import Geometry, XYZGeometry
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

    # Erase old properties and add new event and DFT properties
    geom.global_properties = {}
    geom._update_history('md')
    geom = geom.add_global_properties(res.to_dict())

    # Finish/clean up
    return geom, res

# check lenths block if you have parameters that are global configure at one time


class XTBWrapper(WrapperInterface):
    '''
    Wrapper for xtb functionality.

    Implements :class:`~isicle.interfaces.MDWrapperInterface` to ensure
    required methods are exposed.

    Attributes
    ----------
    temp_dir : str
        Path to temporary directory used for simulation.
    task_map : dict
        Alias mapper for supported molecular dynamic presets. Includes
        "optimize", "crest", "nmr", "protonate", "deprtonate", and "tautomer".
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
        Save molecule
        Parameters
        ----------
        path : str
            Path to input file
            Prefered format is mol or pdb
        method : str 
            GFN forcefield for the optimization
            Default: gff
            Supported forcefields: gfn2, gfn1, gff
        optlevel : str
            Optimization convergence level
            Default : normal
            Supported : crude, sloppy, loose, lax, normal, tight, vtight extreme

        '''

        # Add base command
        s = 'xtb '

        # Add geometry
        #s += os.path.join(self.temp_dir.name,
        #                       '{}.{}'.format(self.geom.basename,
        #                                      self.fmt)

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

        #s += os.path.join(self.temp_dir.name,
        #                       '{}.{}'.format(self.geom.basename,
        #                                      "log"))

        s += '{}.{}'.format(self.geom.basename, "out")
        return s


    def _configure_crest(self, ewin=6, optlevel='Normal', forcefield='gff',
                       protonate=False, deprotonate=False, tautomerize=False,
                       ion=None, charge=None, dryrun=False):

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
        Set up list of xtb jobs

        Parameters
        ----------
        tasks : str
            One task at a time.
        forcefield : str ot list of str
            Forcefield selection. Supply globally or per task.
        ewin : int or list of int
            Energy window for crest calculation. 
        ion : str 
            Ion for protomer calculation.
        optlevel : str or list of str
            Set optimization level. Supply globally or per task. 

        '''

        if type(task) is list:
            raise TypeError('Initiate one xtb or crest job at a time.')
        if type(forcefield) is list:
            raise TypeError('Initiate one forcefield at a time.')
        if type(optlevel) is list:
            raise TypeError('Initiate one opt level at a time.')


        if task is 'optimize':
            config = self._configure_xtb(optlevel=optlevel,
                                         forcefield=forcefield)

        else:
            if task is 'crest':
                p, d, t, i = False, False, False, None

            elif task is 'protonate':
                p, d, t, i = True, False, False, ion

            elif task is 'deprotonate':
                p, d, t, i = False, True, False, ion

            elif task is 'tautomerize':
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
        return self

    def run(self):
        '''
        subprocess.run change to working directory then runs job from command line
        '''
        owd = os.getcwd()
        os.chdir(self.temp_dir.name)
        job = self.config
        os.system(job)
        #with open('job.log', "w") as outfile:
        #    subprocess.run(job, stdout=outfile)

        os.chdir(owd)

    def finish(self, keep_files=False, path=None):
        '''
        Save XTB output files and clean temporary directory

        Parameters 
        ----------
        opt : xtbopt.xyz
        conformer : crest_conformers.xyz and crest.energies, or crest_best.xyz
        protonate : protomers.xyz
        deprotonate : deprotonated.xyz
        tautomers : tautomers.xyz

        '''

        parser = XTBParser()

        parser.load(os.path.join(self.temp_dir.name, self.geom.basename + '.log'))
        result = parser.parse(to_parse=['energy', 'geometry'])

        if keep_files is True:
            import shutil

            if path is None:
                raise ValueError('Must supply `path`.')
            else:
                # TODO: anything else to keep?
                shutil.copy2(os.path.join(self.temp_dir.name,
                                          self.geom.basename + '.out'), path)
                shutil.copy2(os.path.join(self.temp_dir.name,
                                          self.geom.basename + '.log'), path)
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
