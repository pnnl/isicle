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
# TO DO : if optimize is requested along with a crest calculation, need to have the opt 
#         ignore the energy window in cycle
# TO DO : Add in dry run option for adduct to make sure that the adduct of interest is 
#        included in crest


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
    mdw.save_geometry(fmt=kwargs.pop('fmt'))

    # Job type
    mdw.job_type()

    # Run MD simulation
    mdw.run()

    # Create new Geometry with updated structure
    # res['geometry'] will be None or a path to an xyz file.
    geom = geom._update_structure(False, xyz_filename=res['geometry'])

    # Erase old properties and add new event and DFT properties
    geom.global_properties = {}
    geom._update_history('md')
    geom = geom.add_global_properties(res.to_dict())

    # Finish/clean up
    return mdw.finish()

# check lenths block if you have parameters that are global configure at one time


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
        Alias mapper for supported molecular dynamic presets. These include
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

        # All other formats
        self.geom.save(outfile)

        # Path operations
        self.fmt = fmt.lower()
        outfile = os.path.join(self.temp_dir.name,
                               '{}.{}'.format(self.geom.basename,
                                              self.fmt.lower()))


    def load_geometry(self, path):
        # Workaround for xyz input
        # See `isicle.geometry.load_xyz`
        fn, ext = os.path.splitext(path)
        if ext.lower() == '.xyz':
            self.geom = _load_generic_geom(path)
        else:
            self.geom = load(path)

        # Extract filename
        self.geom.basename = os.path.splitext(
            os.path.basename(self.geom.path))[0]

    def set_geometry(self, geom):
        # Assign geometry
        self.geom = geom

        # Extract filename
        self.geom.basename = os.path.splitext(
            os.path.basename(self.geom.path))[0]

    def _job_type_optimize(self, forcefield='gff', optlevel='normal', ewin=False, ion=False, dryrun=False):
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

        d = {'forcefield': forcefield,
             'optlevel': optlevel,
             'basename': self.geom.basename,
             'fmt': self.fmt}

        return ('xtb {basename}.{fmt} --opt {optlevel} --{forcefield} --scratch').format(**d)

    def _job_type_crest(self, forcefield='gff', optlevel='vtight', ewin=6, ion=False, dryrun=False):

        d = {'forcefield': forcefield,
             'optlevel': optlevel,
             'ewin': ewin,
             'basename': self.geom.basename,
             'fmt': self.fmt}

        return ('crest {basename}.{fmt} --ewin {ewin} --optlevel {optlevel} -{forcefield} --scratch').format(**d)

    def _job_type_nmr(self):
        # self.directory = 'xtb_nmr/'
        # TO DO: get protocol for enso/anmr
        raise NotImplementedError

    def _job_type_protonate(self, forcefield='gff', optlevel='normal', ion='H+', ewin=6, dryrun=False):

        if dryrun == True:
            d = {'forcefield': forcefield,
                 'optlevel': optlevel,
                 'ewin': ewin,
                 'ion': ion,
                 'basename': self.geom.basename,
                 'fmt': self.fmt}

            return ('crest {basename}.{fmt} -protonate --swel {ion} --optlevel {optlevel} --ewin {ewin} -{forcefield} --scratch --dryrun').format(**d)

        if dryrun == False:
            d = {'forcefield': forcefield,
                 'optlevel': optlevel,
                 'ewin': ewin,
                 'ion': ion,
                 'basename': self.geom.basename,
                 'fmt': self.fmt}

            return ('crest {basename}.{fmt} -protonate --swel {ion} --optlevel {optlevel} --ewin {ewin} -{forcefield} --scratch').format(**d)

    def _job_type_deprotonate(self, forcefield='gff', optlevel='normal', ewin=6, ion=False, dryrun=False):

        d = {'forcefield': forcefield,
             'optlevel': optlevel,
             'ewin': ewin,
             'basename': self.geom.basename,
             'fmt': self.fmt}

        return ('crest {basename}.{fmt} -deprotonate --optlevel {optlevel} --ewin {ewin} -{forcefield} --scratch').format(**d)

    def _job_type_tautomer(self, forcefield='gff', optlevel='normal', ewin=6, ion=False, dryrun=False):

        d = {'forcefield': forcefield,
             'optlevel': optlevel,
             'ewin': ewin,
             'basename': self.geom.basename,
             'fmt': self.fmt}

        return ('crest {basename}.{fmt} --scratch -{forcefield} --optlevel {optlevel} --ewin {ewin}').format(**d)

    def job_type(self, tasks='optimize', forcefield='gff', ewin=1, ion='H+', optlevel='Normal', dryrun=False):
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
 #       energywindow = safelist(ewin)

        job_list = []

        # Check lengths
        if not ((len(tasks) == len(forcefield)) or (len(forcefield) == 1)):
            raise ValueError('Forcefield must be assigned globally or per'
                             'task.')
        if not ((len(tasks) == len(optlevel)) or (len(optlevel) == 1)):
            raise ValueError('Optimization level must be assigned globally or per'
                             'task.')
#        if not ((len(tasks) == count(ewin)) or (count(ewin) == 1)):
#            raise ValueError('Energy window must be assigned globally or per'
#                             'task.')

# TO DO: cycle so ewin and ion only goes to ''
        for task, f, o in zip(tasks, cycle(forcefield), cycle(optlevel)):
            job_list.append(self.task_map[task](forcefield=f,
                                                optlevel=o,
                                                ewin=ewin,
                                                ion=ion,
                                                dryrun=dryrun))

        self.job_list = job_list

        return self.job_list

    def save_job_type(self):
        '''
        Write 

        '''

        # Write to file
        with open(os.path.join(self.temp_dir.name,
                               self.geom.basename + '.nw'), 'w') as f:
            f.write(self.config)

    def run(self):
        '''
        subprocess.run change to working directory then runs job from command line

        iterate jobs from job_type. have new temp_dir for each job?

        '''

        os.chdir(self.temp_dir.name)

        for i in self.job_list:
            print(i)
            subprocess.call(i, shell=True)

    def finish(self, keep_files=False, path=None):
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

        if keep_files is True:
            import shutil
            import glob

            if path is None:
                raise ValueError('Must supply `path`.')
            else:
                # TODO: anything else to keep?
                if 'xtbopt.xyz' in glob.glob(self.temp_dir):
                    shutil.copy2(os.path.join(self.temp_dir.name,
                                              'xtbopt.xyz'), path)
                if 'crest*' in glob.glob(self.temp_dir):
                    shutil.copy2(os.path.join(self.temp_dir.name,
                                              'crest.energies'), path)
                    shutil.copy2(os.path.join(self.temp_dir.name,
                                              'crest_conformers.xyz'), path)
                    shutil.copy2(os.path.join(self.temp_dir.name,
                                              'crest_best.xyz'), path)
                if 'prot*' in glob.glob(self.temp_dir):
                    shutil.copy2(os.path.join(self.temp_dir.name,
                                              'protomers.xyz'), path)
                if 'deprot*' in glob.glob(self.temp_dir):
                    shutil.copy2(os.path.join(self.temp_dir.name,
                                              'deprotonated.xyz'), path)
                if 'taut*' in glob.glob(self.temp_dir):
                    shutil.copy2(os.path.join(self.temp_dir.name,
                                              'tautomers.xyz'), path)

        # Remove temporary files
        self.temp_dir.cleanup()


def amber():
    # Don't work on this one yet
    # TODO: define input, arguments, and output
    # Should return instance (or list) of MDOptimizedGeometry
    raise NotImplementedError

# TODO: add other MD methods as needed
