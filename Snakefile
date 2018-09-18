from os.path import *
from resources.utils import cycles, frames

# snakemake configuration
configfile: 'config.yaml'
localrules: all
# include: 'impact.snakefile'
include: 'molecular_dynamics.snakefile'


# a pseudo-rule that collects the target files

# # impact
# rule all:
#     input:
#         join(config['path'], 'output', 'impact_results.tsv')


# through md
IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))

rule all:
    input:
        expand(join(config['path'], 'output', 'selected', 'xyz', '{id}_{adduct}_{cycle}_s.xyz'),
               id=IDS, adduct=config['adducts'], cycle=cycles()[1:])
