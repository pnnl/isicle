from os.path import *
from resources.nwchem.parseOutputs import XYZtoMFJ

# snakemake configuration
configfile: 'isicle/config.yaml'
include: 'molecular_dynamics.snakefile'


rule copyOver:
    input:
        join(config['path'], 'output', 'md', 'downselected', '{id}_Ne_{cycle}_{selected}.xyz')
    output:
        join(config['path'], 'output', 'shielding', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.xyz')
    # group:
    #     'shielding'
    shell:
        'cp {input} {output}'

# create .nw files based on template (resources/nwchem/template.nw)
rule createNW:
    input:
        rules.copyOver.output
    output:
        join(config['path'], 'output', 'shielding', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.nw')
    # group:
    #     'shielding'
    shell:
        'python isicle/resources/nwchem/generateNW.py {input} --template {config[nwchem][shielding_template]}'

# run NWChem
rule NWChem:
    input:
        xyz = rules.createNW.input,
        nw = rules.createNW.output
    output:
        join(config['path'], 'output', 'shielding', '{id}', 'cycle_{cycle}_{selected}', '{id}_{cycle}_{selected}.out')
    benchmark:
        join(config['path'], 'output', 'shielding', 'benchmarks', '{id}_{cycle}_{selected}.nwchem.benchmark')
    # group:
    #     'shielding'
    shell:
        '{config[nwchem][runscript]} {input.nw}'
