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
        # expand(join(config['path'], 'output', 'chemical_shifts', '{id}.tsv'), id=IDS)
        expand(join(config['path'], 'output', 'shielding', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.nw'),
               id=IDS, cycle=cycles(config['amber']['cycles']), selected=['s', 'd1', 'd2'])
