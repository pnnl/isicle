from os.path import *
from isicle.utils import cycles
import shutil
from pkg_resources import resource_filename

# snakemake configuration
include: 'shielding.snakefile'

SMI, = glob_wildcards(join(config['path'], 'input', '{id}.smi'))
INCHI, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))
IDS = SMI + INCHI

# copy reference molecule
if config['nwchem']['reference'] in ['TMS', 'DSS']:
    shutil.copy2(resource_filename('isicle', join('resources', 'nwchem', config['nwchem']['reference'] + '.smi')),
                 join(config['path'], 'input'))
    IDS.append(config['nwchem']['reference'])

# check if supplied reference is in input
elif config['nwchem']['reference'] not in IDS:
    raise Exception('Select TMS or DSS as reference, or ensure %s is in the input folder.' % config['nwchem']['reference'])


rule all:
    input:
        expand(join(config['path'], 'output', 'shifts', '{id}.tsv'), id=IDS)
