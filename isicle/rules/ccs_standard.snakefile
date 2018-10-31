from os.path import *

# snakemake configuration
include: 'mobility.snakefile'

SMI, = glob_wildcards(join(config['path'], 'input', '{id}.smi'))
INCHI, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))
IDS = SMI + INCHI


rule all:
    input:
        expand(join(config['path'], 'output', 'mobility', 'mobcal', 'calibrated_ccs', '{id}_{adduct}.tsv'),
               id=IDS, adduct=config['adducts'])
