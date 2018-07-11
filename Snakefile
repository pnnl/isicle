from os.path import *

# snakemake configuration
configfile: 'config.yaml'
localrules: all
include: 'impact.snakefile'

# infer wildcards from inputs, select adducts
IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))

# a pseudo-rule that collects the target files
rule all:
    input:
        # expand(join(config['path'], 'output', '5_impact', '{id}_{adduct}.txt'),
        #        id=IDS, adduct=config['adducts'])
        join(config['path'], 'output', 'impact_results.tsv')
