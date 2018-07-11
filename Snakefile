from os.path import *

# snakemake configuration
configfile: 'config.yaml'
localrules: all
include: 'impact.snakefile'

# a pseudo-rule that collects the target files
rule all:
    input:
        join(config['path'], 'output', 'impact_results.tsv')
