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

# TO DO: Have finish() save xyz contents to a callable object in memory
# TO DO: Set up generalized temp directory to run xtb jobs in
# TO DO: Finish job_type(), possibly combine with run?

def _program_selector(program):
    program_map = {'xtb': XTBWrapper}

    if program.lower() in program_map.keys():
        return program_map[program.lower()]()
    else:
        raise ValueError('{} not a supported molecular dyanmics program.'.format(program))

def md(self, path, program='xtb', template=None, **kwargs):
    '''
    Optimize geometry, either XYZ or PDB, using stated program.
    Additional inputs can be grid size, optimization criteria level,
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

class XTBWrapper(Geometry, MDWrapperInterface):
    def __init__(self,):
        self.directory = None
        self.commandline = None
        self.temp_dir = tempfile.TemporaryDirectory()
        self.task_map = {'optimize': self._job_type_optimize,
                         'crest': self._job_type_crest,
                         'nmr': self._job_type_nmr,
                         'protonate': self._job_type_protonate,
                         'deprotonate': self._job_type_deprotonate,
                         'tautomer': self._job_type_tautomer}

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
        directory : str
            Directory where job will be run
        method : str 
            GFN forcefield for the optimization
            Default: gff
            Supported forcefields: gfn2, gfn1, gff
        optlevel : str
            Optimization convergence level
            Default : normal
            Supported : crude, sloppy, loose, lax, normal, tight, vtight extreme
            Format to save this molecule in. If None, determined from given
            path's extension. If .pkl. pickles this full object.
            Default: None.
            Supported formats: .smi (SMILES), .inchi (InChI), .mol, .xyz,
            .pdb, .pkl.
        Returns
        -------
        xyz file
        '''

        self.directory = self.filename
        self.commandline = 'xtb {input} --opt {optlevel} --{method}'

    def _job_type_crest(self, input, forcefield='gff', optlevel='vtight', ewin=1):
        '''
        Runs Conformer Rotamer Ensemble Sampling Tool
        Parameters
        ----------
        path : str
            Path to input file
        directory : str
            Directory where job will be run
        method : str 
            GFN forcefield for the optimization
            Default: gff
            Supported forcefields: gfn2, gfn1, gff
        optlevel : str
            Optimization convergence level
            Default : vtight
            Supported : crude, sloppy, loose, lax, normal, tight, vtight extreme
        ewin : int
            Set the energy threshold to int kcal/mol
            Default : 6 kcal/mol
        Returns
        -------
        xyz file(s)  
        '''

        self.directory = 'crest/' + self.geom_basename
        self.commandline = ''

        return 

    def _job_type_nmr(self, input):
        # self.directory = 'xtb_nmr/'
        # TO DO: get protocol for enso/anmr
        raise NotImplementedError

    def _job_type_protonate(self, input, ion='H+'):
        return ['protonate/', 'crest {input} -protonate --swel {ion}']

    def _job_type_deprotonate(self, input):
        return ['deprotonate/', 'crest {input} -deprotonate']

    def _job_type_tautomer(self, input):
        self.directory = 'tautomer/'
        self.commandline = 'crest {input}'

    def job_type(self, input, directory = self.directory, tasks='optimize', forcefield='gfn2'):
        '''
        Returns the directory and command line(s) needed to run md job 
        '''
        
        # Cast to list safely
        tasks = safelist(tasks)
        forcefield = safelist(forcefield)


        subprocess.run('mkdir {directory}')
        job_list =[]

        for task, f, e in zip(tasks, cycle(forcefield)):
            # TODO: finish this
            config += self.task_map[task](functional=f, ewin=)


    def run(self, opt=False, crest=False, enso=False):
        '''
        subprocess.run change to working directory then runs job from command line
        loop through an array 
        '''
        subprocess.run()



    def finish(self, keep_files=True, path=None):
        '''
        Saves
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
                shutil.copy2(os.path.join(self.temp_dir.name, self.geom.basename + 'xtbopt.xyz'), path)
                shutil.copy2(os.path.join(self.temp_dir.name, self.geom.basename + 'crest.energies'), path)
                shutil.copy2(os.path.join(self.temp_dir.name, self.geom.basename + 'crest_conformers.xyz'), path)
                shutil.copy2(os.path.join(self.temp_dir.name, self.geom.basename + 'crest_best.xyz'), path)
                shutil.copy2(os.path.join(self.temp_dir.name, self.geom.basename + 'protomers.xyz'), path)
                shutil.copy2(os.path.join(self.temp_dir.name, self.geom.basename + 'deprotonated.xyz'), path)

                geoms = glob.glob(os.path.join(self.temp_dir.name, '*.{}'.format(self.fmt)))
                [shutil.copy2(x, path) for x in geoms]

        # Remove temporary files
        self.temp_dir.cleanup()


class AMBERWrapper(Geometry, MDWrapperInterface)


    def amber():
        # Don't work on this one yet
        # TODO: define input, arguments, and output
        # Should return instance (or list) of MDOptimizedGeometry
        raise NotImplementedError

# TODO: add other MD methods as needed
