from os.path import *
from glob import glob

# snakemake configuration
configfile: 'config.yaml'
localrules: all
include: 'impact.snakefile'

# # infer wildcards from inputs, select adducts
IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))

# a pseudo-rule that collects the target files
rule all:
    input:
        expand(join(config['path'], 'output', '5_impact', '{id}_{adduct}.txt'),
               id=IDS, adduct=config['adducts'])
