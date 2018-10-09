from os.path import *

# snakemake configuration
include: 'isicle/rules/mobility.snakefile'

# end to end
IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))

rule all:
    input:
        expand(join(config['path'], 'output', 'mobility', 'mobcal', 'boltzmann_ccs', '{id}_{adduct}.tsv'), id=IDS, adduct=config['adducts'])
