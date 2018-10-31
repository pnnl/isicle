from os.path import *

# snakemake configuration
include: 'mobility_alt.snakefile'

SMI, = glob_wildcards(join(config['path'], 'input', '{id}.smi'))
INCHI, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))
IDS = SMI + INCHI


rule all:
    input:
        expand(join(config['path'], 'output', 'mobility', 'impact', 'ccs', '{id}_{adduct}.He.ccs'),
               id=IDS, adduct=config['adducts']),
        expand(join(config['path'], 'output', 'mobility', 'impact', 'ccs', '{id}_{adduct}.N2.ccs'),
               id=IDS, adduct=config['adducts'])
