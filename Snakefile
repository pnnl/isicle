from os.path import *
from resources.utils import *

# snakemake configuration
configfile: 'config.yaml'

# include: 'impact.snakefile'
# include: 'molecular_dynamics.snakefile'
# include: 'dft.snakefile'
include: 'mobility.snakefile'


# # impact
# IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))

# rule all:
#     input:
#         expand(join(config['path'], 'output', '5_impact', '{id}_{adduct}.N2.ccs'), id=IDS, adduct=config['adducts'])


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

# # end to end
# IDS, = glob_wildcards(join(config['path'], 'input', '{id}.inchi'))

rule all:
    input:
        expand(join(config['path'], 'output', 'boltzmann_ccs', '{id}_{adduct}.tsv'), id=IDS, adduct=config['adducts'])
