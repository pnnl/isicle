from os.path import *

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


def cycles():
    return ['%03d' % x for x in range(config['cycles'] + 1)]


def frames():
    return ['%03d' % x for x in range(config['nframes'])]


rule all:
    input:
        expand(join(config['path'], 'output', 'selected', 'xyz', '{id}_{adduct}_{cycle}_s.xyz'),
               id=IDS, adduct=config['adducts'], cycle=cycles()[1:])
