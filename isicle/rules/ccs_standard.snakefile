from os.path import *

# snakemake configuration
include: 'mobility.snakefile'

# default to SMILES as input
IDS, = glob_wildcards(join(config['path'], 'input', '{id}.smi'))

# if none found, check for InChIs
if len(IDS) == 0:
    IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))


rule all:
    input:
        expand(join(config['path'], 'output', 'mobility', 'mobcal', 'boltzmann_ccs', '{id}_{adduct}.tsv'),
               id=IDS, adduct=config['adducts'])
