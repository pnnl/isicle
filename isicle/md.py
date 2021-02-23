from isicle.interfaces import MDWrapperInterface
from isicle import geometry
import subprocess
import tempfile
import os
from isicle.utils import safelist
from itertools import combinations, cycle

'''
Files resulting from an xtb job always run in the same directory that the command is 
issued in, no matter where the input is. Can direct the .log file, but no other files.
'''

# TO DO : Add implicit solvation, frequency
# TO DO : if optimize is requested along with a crest calculation, need to have the opt ignore the energy window in cycle
# TO DO : Add in dry run option for adduct to make sure that the adduct of interest is included in crest

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
        raise ValueError('{} not a supported molecular dyanmics program.'.format(program))

def md(self, path, program='xtb', **kwargs):
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

    # Load geometry
    mdw.load_geometry(path)

    # Save geometry
    mdw.save_geometry(path, fmt=kwargs.pop('fmt'))

    # Job type
    mdw.job_type()

    # Run MD simulation
    mdw.run()

    # Finish/clean up
    return mdw.finish()

#check lenths block if you have parameters that are global configure at one time

class XTBWrapper(MDWrapperInterface):
    '''
    Wrapper for xtb functionality.

    Implements :class:`~isicle.interfaces.MDWrapperInterface` to ensure
    required methods are exposed.

    Attributes
    ----------
    temp_dir : str
        Path to temporary directory used for simulation.
    task_map : dict
        Alias mapper for supported molecular dynamic presets. Thses include
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
        self.task_map = {'optimize': self._job_type_optimize,
                         'crest': self._job_type_crest,
                         'nmr': self._job_type_nmr,
                         'protonate': self._job_type_protonate,
                         'deprotonate': self._job_type_deprotonate,
                         'tautomer': self._job_type_tautomer}

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

        # Extract filename
        self.geom.basename = os.path.splitext(os.path.basename(self.geom.path))[0]

    def save_geometry(self, fmt='xyz'):
        '''
        Save internal :obj:`~isicle.geometry.Geometry` representation to file.

        Parameters
        ----------
        fmt : str
            Filetype used by xtb. Must be "xyz", "smi", ".inchi", ".mol", ".xyz",
            ".pdb", ".pkl".

        '''

        # Path operations
        self.fmt = fmt.lower()
        outfile = os.path.join(self.temp_dir.name,
                               '{}.{}'.format(self.geom.basename,
                                              self.fmt.lower()))

        # Workaround for xyz input
        # See `isicle.geometry.load_xyz`
        if self.geom.filetype == '.xyz':
            if self.fmt != 'xyz':
                raise TypeError('Input .xyz files cannot be converted.')

            with open(outfile, 'w') as f:
                f.write('\n'.join(self.geom.contents))
            return

        # All other formats
        self.geom.save(outfile)

    def load_geometry(self, path):
        # Workaround for xyz input
        # See `isicle.geometry.load_xyz`
        fn, ext = os.path.splitext(path)
        if ext.lower() == '.xyz':
            self.geom = _load_generic_geom(path)
        else:
            self.geom = load(path)

        # Extract filename
        self.geom.basename = os.path.splitext(os.path.basename(self.geom.path))[0]

    def set_geometry(self, geom):
        # Assign geometry
        self.geom = geom

        # Extract filename
        self.geom.basename = os.path.splitext(os.path.basename(self.geom.path))[0]

    def _job_type_optimize(self, input, forcefield='gff', optlevel='normal'):

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
    
        Returns
        -------
        xyz file
        '''

        return 'xtb {self.geom} --opt {optlevel} --{method} --scratch'

    def _job_type_crest(self, input, forcefield='gff', optlevel='vtight', ewin=6):

        return 'crest {self.geom} --ewin {ewin} --optlevel {optlevel} -{forcefield} --scratch'

    def _job_type_nmr(self, input):
        # self.directory = 'xtb_nmr/'
        # TO DO: get protocol for enso/anmr
        raise NotImplementedError

    def _job_type_protonate(self, input, ion='H+', ewin=6, dryrun=False):

        return 'crest {self.geom} -protonate --swel {ion} --ewin {ewin} --scratch'

    def _job_type_deprotonate(self, input, ewin=6):

        return 'deprotonate/', 'crest {self.geom} -deprotonate --ewin {ewin} --scratch'

    def _job_type_tautomer(self, input, ewin=6):

        return 'crest {self.geom} --scratch, --ewin {ewin}'


    def job_type(self, tasks='optimize', forcefield='gff', ewin=1, ion='H+', optlevel='Normal',dryrun=False):
        '''
        Set up list of xtb jobs

        Parameters
        ----------
        tasks : str or list of str
            Tasks text.
        forcefield : str ot lisst of str
            Forcefield selection. Supply globally or per task.
        ewin : int or list of int
            Energy window for crest calculation. 
        ion : str 
            Ion for protomer calculation.
        optlevel : str or list of str
            Set optimization level. Supply globally or per task. 

        '''
        
        # Cast to list safely
        tasks = safelist(tasks)
        optlevel = safelist(optlevel)
        forcefield = safelist(forcefield)
        energywindow = safelist(energywindow)

        job_list =[]


        # Check lengths
        if not ((len(tasks) == len(forcefield)) or (len(functional) == 1)):
            raise ValueError('Forcefield must be assigned globally or per'
                             'task.')
        if not ((len(tasks) == len(optlevel)) or (len(optlevel) == 1)):
            raise ValueError('Optimization level must be assigned globally or per'
                             'task.')
        if 'optimize' in tasks:   
            if not ((len(tasks) == len(ewin)) or (len(ewin) == 1)):
                raise ValueError('Energy window must be assigned globally or per'
                                 'task.')

        for task, f, o, e in zip(tasks, cycle(forcefield), cycle(optlevel), cycle(ewin)):
            job_list.append(self.task_map[task](forcefield=f, 
                                                ewin=ewin,
                                                optlevel=o,
                                                ion=ion))

        self.job_list = job_list

        return self.job_list


    def run(self, keep_files=True, path=None):
        '''
        subprocess.run change to working directory then runs job from command line

        iterate jobs from job_type. have new temp_dir for each job?

        '''

        os.chdir(self.temp_dir.name)

        for i in self.job_list:
            subprocess.run(i)


    def finish(self, keep_files=True, path=None):
        '''
        Save XTB output files and clean temporary direc

        Parameters 
        -------
        opt: xtbopt.xyz
        conformer: crest_conformers.xyz and crest.energies, or crest_best.xyz
        enso/anmr: anmr.dat or newanmr.dat
        protonate: protomers.xyz
        deprotonate: deprotonated.xyz
        tautomers: tautomers.xyz

        '''

        #TO DO : Option if the file isn't abailale? I.e. you didn't run that job or the job failed.

        if keep_files is True:
            import shutil
            import glob

            if path is None:
                raise ValueError('Must supply `path`.') 
            else:
                # TODO: anything else to keep?
                shutil.copy2(os.path.join(self.temp_dir.name, 
                                          'xtbopt.xyz'), path)
                shutil.copy2(os.path.join(self.temp_dir.name, 
                                          'crest.energies'), path)
                shutil.copy2(os.path.join(self.temp_dir.name, 
                                          'crest_conformers.xyz'), path)
                shutil.copy2(os.path.join(self.temp_dir.name, 
                                          'crest_best.xyz'), path)
                shutil.copy2(os.path.join(self.temp_dir.name, 
                                          'protomers.xyz'), path)
                shutil.copy2(os.path.join(self.temp_dir.name, 
                                          'deprotonated.xyz'), path)

                geoms = glob.glob(os.path.join(self.temp_dir.name, '*.{}'.format(self.fmt)))
                [shutil.copy2(x, path) for x in geoms]

        # Remove temporary files
        self.temp_dir.cleanup()

def amber():
    # Don't work on this one yet
    # TODO: define input, arguments, and output
    # Should return instance (or list) of MDOptimizedGeometry
    raise NotImplementedError

# TODO: add other MD methods as needed
