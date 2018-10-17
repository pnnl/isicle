from os.path import *
from isicle.utils import cycles

# snakemake configuration
include: 'shielding.snakefile'

# default to SMILES as input
IDS, = glob_wildcards(join(config['path'], 'input', '{id}.smi'))

# if none found, check for InChIs
if len(IDS) == 0:
    IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))


rule all:
    input:
        expand(join(config['path'], 'output', 'shielding', 'boltzmann_shielding', '{id}.tsv'), id=IDS)
