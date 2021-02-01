from isicle import geometry
import os
import 

'''
Files resulting from an xtb job always run in the same directory that the command is 
issued in, no matter where the input is. Can direct the .log file, but no other files.
'''

def xtb_optimize(path, directory='xtb_output/', method='gff', optlevel='normal'):
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

    # TODO: define output
    # Load file, can be mol or xyz file if xyz we need a .chrg file

    # Save file to working directory
    input = geometry.save(path, directory, fmt='mol')

    cmd = 'xtb {input} --opt {optlevel} --{method} '

    # Change to work directory and submit job
    os.chdir(directory)
    os.system(cmd)

    opt_structure = directory+'/xtbopt.xyz'

    return opt_structure


def crest(path, directory='crest_output/', method='gff', optlevel='vtight', ewin=1):
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

    # TODO: Add option to use crest.energies to downselect conformers, otherwise return all conformers.

    # Save file to working directory
    input = geometry.save(path, directory, fmt='mol')
    cmd = 'crest {input} --opt {optlevel} -{method} -ewin {ewin}'
 
    # Change to work directory and submit job
    os.chdir(directory)
    os.system(cmd)

    # Split crest_conformers.xyz into separate xyz files
    fname = 'split/' + (os.path.splitext(path)[0]).split('/')[-1] + '_conf_'
    split = 'obabel crest_conformers.xyz -O {fname} -m'
    os.system(split)

    # Option to downselect
    
    raise NotImplementedError


def enso(path, directory, method='gff', optlevel='vtight'):
    # TODO: define input, arguments, and output
    # Should return instance (or list) of ENSO calculated chemical shifts
    cmd = 'crest {input} -{method} -g chcl3 -T 4 -nmr > crest.out'
    raise NotImplementedError


def amber():
    # Don't work on this one yet
    # TODO: define input, arguments, and output
    # Should return instance (or list) of MDOptimizedGeometry
    raise NotImplementedError

# TODO: add other MD methods as needed
