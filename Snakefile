from os.path import *
from isicle.resources.utils import *

# snakemake configuration
configfile: 'isicle/config.yaml'
include: 'isicle/mobility.snakefile'

# end to end
IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))

rule all:
    input:
        expand(join(config['path'], 'output', 'boltzmann_ccs', '{id}_{adduct}.tsv'), id=IDS, adduct=config['adducts'])
