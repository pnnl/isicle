from os.path import *

# snakemake configuration
include: 'mobility_alt.snakefile'

# default to SMILES as input
IDS, = glob_wildcards(join(config['path'], 'input', '{id}.smi'))

# if none found, check for InChIs
if len(IDS) == 0:
    IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))


rule all:
    input:
        expand(join(config['path'], 'output', 'mobility', 'impact', 'ccs', '{id}_{adduct}.He.ccs'),
               id=IDS, adduct=config['adducts']),
        expand(join(config['path'], 'output', 'mobility', 'impact', 'ccs', '{id}_{adduct}.N2.ccs'),
               id=IDS, adduct=config['adducts'])
