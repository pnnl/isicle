from os.path import *
from isicle.utils import cycles
import shutil
from pkg_resources import resource_filename

# snakemake configuration
include: 'shielding.snakefile'

# default to SMILES as input
IDS, = glob_wildcards(join(config['path'], 'input', '{id}.smi'))

# if none found, check for InChIs
if len(IDS) == 0:
    IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))

# copy reference molecule
if config['nwchem']['reference'] in ['TMS']:
    if not exists(join(config['path'], 'output', 'adducts', 'geometry_Ne', config['nwchem']['reference'] + '_Ne.mol2')):
        shutil.copy2(resource_filename('isicle', join('resources', 'nwchem', config['nwchem']['reference'] + '.mol2')),
                     join(config['path'], 'output', 'adducts', 'geometry_Ne', config['nwchem']['reference'] + '_Ne.mol2'))
    IDS.append(config['nwchem']['reference'])


rule all:
    input:
        expand(join(config['path'], 'output', 'shifts', '{id}.tsv'), id=IDS)
