from os.path import *

# snakemake configuration
include: 'mobility.snakefile'

SMI, = glob_wildcards(abspath(join('input', '{id}.smi')))
INCHI, = glob_wildcards(abspath(join('input', '{id}.inchi')))
IDS = SMI + INCHI

IDS.sort()

if 'stop' in config:
    IDS = IDS[config['start']:config['stop']]


rule all:
    input:
        expand(abspath(join('output', 'mobility', 'mobcal', 'calibrated_ccs', '{id}_{adduct}.tsv')),
               id=IDS, adduct=config['adducts'])
