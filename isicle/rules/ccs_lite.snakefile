from os.path import *

# snakemake configuration
include: 'mobility_alt.snakefile'

SMI, = glob_wildcards(join('input', '{id}.smi'))
INCHI, = glob_wildcards(join('input', '{id}.inchi'))
IDS = SMI + INCHI

IDS.sort()
IDS = IDS[config['start']:config['stop']]


rule all:
    input:
        expand(join('output', 'mobility', 'impact', 'ccs', '{id}_{adduct}.He.ccs'),
               id=IDS, adduct=config['adducts']),
        expand(join('output', 'mobility', 'impact', 'ccs', '{id}_{adduct}.N2.ccs'),
               id=IDS, adduct=config['adducts'])
