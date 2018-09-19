from os.path import *
from resources.utils import *

# snakemake configuration
configfile: 'config.yaml'
localrules: all
# include: 'impact.snakefile'
# include: 'molecular_dynamics.snakefile'
# include: 'dft.snakefile'
include: 'mobility.snakefile'


# a pseudo-rule that collects the target files

# # impact
# rule all:
#     input:
#         join(config['path'], 'output', 'impact_results.tsv')


# # through md
# IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))

# rule all:
#     input:
#         expand(join(config['path'], 'output', 'selected', 'xyz', '{id}_{adduct}_{cycle}_{selected}.xyz'),
#                id=IDS, adduct=config['adducts'], cycle=cycles(config['cycles'])[1:], selected=['s', 'd1', 'd2'])

# # through dft
# IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))

# rule all:
#     input:
#         expand(join(config['path'], 'output', 'mobcal', '{id}_{adduct}_{cycle}_{selected}_geom+charge.mfj'),
#                id=IDS, adduct=config['adducts'], cycle=cycles(config['cycles'])[1:], selected=['s', 'd1', 'd2'])

# # through mobcal
# IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))

# rule all:
#     input:
#         expand(join(config['path'], 'output', 'mobcal', '{id}_{adduct}_{cycle}_{selected}_geom+charge.out'),
#                id=IDS, adduct=config['adducts'], cycle=cycles(config['cycles'])[1:], selected=['s', 'd1', 'd2'])

# end to end
rule all:
    input:
        join(config['path'], 'output', 'ccs_result.tsv')
